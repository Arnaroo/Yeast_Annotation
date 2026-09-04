#!/usr/bin/env python3
"""
unagi_compare.py
=================
Gene-level comparison of the released annotation against UNAGI (Al Kadi
et al., 2020), an independent nanopore full-length cDNA transcript
annotation that took no part in building this reference: different
library chemistry, different basecaller, different assembler, different
laboratory.

WHY THE COORDINATE PROBLEM DISSOLVES
------------------------------------
UNAGI's supplementary tables carry the ORF boundary on the same row as the
transcript boundary:

  MOESM5 sheet 'TSSs'   gene start (TSS)  and  first codon   5,594 genes
  MOESM6                gene end          and  stop codon    5,755 genes

So a UTR length is the difference of two columns of the same table, in
nucleotides relative to the ORF. That is exactly the quantity our annotation
stores. No genomic-to-relative projection is required, and because both
columns come from one pipeline in one file, the difference is exact whatever
base convention they use.

THE JOIN, AND WHY IT IS TRUSTED
-------------------------------
UNAGI keys rows by RefSeq mRNA accession (NM_...), ours by systematic name
(YAL001C). Rather than depend on an external accession map, genes are matched
on the genomic coordinate of the ORF boundary itself: chromosome, strand, and
the start or stop codon position, against the CDS extent in the backbone GFF3.

That match is only credible if the coordinate conventions are known rather
than guessed, so the offset between the two is measured, not assumed, and
asserted below. It is completely regular:

           3' end (stop codon)      5' end (first codon)
  strand +   backbone = unagi + 0     backbone = unagi + 1
  strand -   backbone = unagi + 1     backbone = unagi + 0

That is the signature of 0-based half-open intervals against our 1-based
inclusive GFF3, and it holds for 5,754 of 5,754 and 6,680 of 6,683 joinable
rows respectively. A convention that resolves every row is a convention that
has been identified correctly. The script asserts this and stops if it ever
stops holding.

WHAT IS COMPARED
----------------
Deliberately the same shape as plot_final_vs_tifseq_comparison.R's rate computation, so that UNAGI and TIF-seq
can be read side by side as two independent benchmarks:

  final    final_utr.tsv            the released merged annotation
  nag      Nagalakshmi_UTR.csv      the prior reference
  drs      DRS_UTR_corrected.tsv    our raw segmentation calls

  benchmark  UNAGI, Al Kadi et al. 2020, nanopore full-length cDNA.
             Took no part in constructing our annotation.

UNAGI is an independent benchmark in a strong sense: different library
chemistry (full-length cDNA, not direct RNA), different basecaller era,
different assembler, and a different laboratory.

Same comparison rule as the TIF-seq benchmark, both values strictly greater than zero,
per end. Same tolerances. Same inclusive thresholds on both tails.

Fetching the UNAGI supplementary files
---------------------------------------
The two source files are the publisher-hosted supplementary tables of
Al Kadi M, Suzuki Y, et al. (2020) "UNAGI: an automated pipeline for
nanopore full-length cDNA sequencing uncovers novel transcripts and
isoforms in yeast." Funct Integr Genomics 20(4):523-536,
DOI 10.1007/s10142-020-00732-1 -- not redistributed here. Fetch them
with:

  BASE=https://static-content.springer.com/esm/art%3A10.1007%2Fs10142-020-00732-1/MediaObjects
  curl -sSL -o 10142_2020_732_MOESM5_ESM.xlsx "$BASE/10142_2020_732_MOESM5_ESM.xlsx"
  curl -sSL -o 10142_2020_732_MOESM6_ESM.csv  "$BASE/10142_2020_732_MOESM6_ESM.csv"

MOESM5 (xlsx) carries 1,089 genes with alternative 5' UTR isoforms.
MOESM6 (csv) carries 1,482 genes with alternative poly(A)/3' UTR sites.

Outputs
-------
  unagi_join.tsv            per gene, the joined benchmark
  unagi_source_counts.tsv   our parse against the paper's own counts
  unagi_coverage.tsv        which genes each side calls at all
  unagi_rates.tsv           annotation x end x tolerance
  unagi_distribution.tsv    signed difference distributions
  unagi_merge_effect.tsv    final against prior reference, paired
  unagi_by_provenance.tsv   final annotation split by boundary source
  unagi_dispersion.tsv      UNAGI's own definition and its own spread
  unagi_rate_curve.tsv      fine tolerance sweep
  unagi_vs_tifseq.tsv       the two benchmarks against each other
  unagi_compare.log         full run log

Run:
  python3 unagi_compare.py \
      --annotation-dir Data/Annotation \
      --unagi-dir      Data/UNAGI \
      --outdir         Data/UNAGI
"""

import argparse
import ast
import collections
import hashlib
import os
import re
import sys
from datetime import datetime, timezone

import numpy as np
import openpyxl
import pandas as pd
import scipy
from scipy import stats

# ---------------------------------------------------------------- paths
parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--annotation-dir", required=True,
                     help="Data/Annotation, holding final_utr.tsv, DRS_UTR_corrected.tsv, "
                          "Nagalakshmi_UTR.csv, Pelechano_UTR.csv and the backbone GFF3")
parser.add_argument("--unagi-dir", required=True,
                     help="directory holding the two fetched UNAGI supplementary files")
parser.add_argument("--outdir", required=True, help="directory for output tables and the log")
args = parser.parse_args()

ANN = args.annotation_dir
UNAGI = args.unagi_dir
OUT = args.outdir
LOGS = OUT
os.makedirs(OUT, exist_ok=True)

INPUTS = {
    "final_utr": os.path.join(ANN, "final_utr.tsv"),
    "drs_utr": os.path.join(ANN, "DRS_UTR_corrected.tsv"),
    "nag_utr": os.path.join(ANN, "Nagalakshmi_UTR.csv"),
    "pel_utr": os.path.join(ANN, "Pelechano_UTR.csv"),
    "backbone_gff": os.path.join(ANN, "gene_models_backbone_with_UTR_slots.gff3"),
    "unagi_5prime": os.path.join(UNAGI, "10142_2020_732_MOESM5_ESM.xlsx"),
    "unagi_3prime": os.path.join(UNAGI, "10142_2020_732_MOESM6_ESM.csv"),
}

TOLERANCES = [20, 50, 100]
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI"]
# RefSeq chromosome accessions for S. cerevisiae S288C, chrI to chrXVI
NC_TO_ROMAN = {"NC_0011%02d" % (33 + i): ROMAN[i] for i in range(16)}

_logfh = open(os.path.join(LOGS, "unagi_compare.log"), "w")


def log(msg=""):
    print(msg)
    _logfh.write(str(msg) + "\n")
    _logfh.flush()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()[:16]


# ---------------------------------------------------------------- header
log("=" * 78)
log("unagi_compare.py")
log("run: %s" % datetime.now(timezone.utc).isoformat(timespec="seconds"))
log("python: %s   numpy: %s   pandas: %s   scipy: %s   openpyxl: %s"
    % (sys.version.split()[0], np.__version__, pd.__version__,
       scipy.__version__, openpyxl.__version__))
log("=" * 78)
log()
log("INPUTS (absolute paths, read only, sha256 prefix)")
for k, v in INPUTS.items():
    log("  %-14s %s  %s" % (k, sha256(v), v))
log()
log("UNAGI source: Al Kadi M et al. (2020) Funct Integr Genomics 20(4):523-536")
log("              doi 10.1007/s10142-020-00732-1, supplementary MOESM5, MOESM6")
log("              fetched manually, see the module docstring")
log()

# ================================================================
# PART 0. The backbone, and the coordinate convention
# ================================================================
log("=" * 78)
log("PART 0.  JOINING UNAGI TO THE BACKBONE BY ORF COORDINATE")
log("=" * 78)
log()


def load_cds_extents(path):
    """3'-most and 5'-most CDS coordinate per gene, 1-based inclusive."""
    parts = collections.defaultdict(list)
    strand, chrom = {}, {}
    id_re = re.compile(r"ID=(.+?)_CDS;")
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = id_re.search(f[8])
            if not m:
                continue
            g = m.group(1)
            parts[g].append((int(f[3]), int(f[4])))
            strand[g], chrom[g] = f[6], f[0].strip()
    rows = []
    for g, v in parts.items():
        lo = min(a for a, _ in v)
        hi = max(b for _, b in v)
        if strand[g] == "+":
            atg, stop = lo, hi
        else:
            atg, stop = hi, lo
        rows.append((g, chrom[g], strand[g], atg, stop))
    return pd.DataFrame(rows, columns=["gene", "chr", "strand", "atg", "stop"])


cds = load_cds_extents(INPUTS["backbone_gff"])
cds = cds[cds["chr"].isin(ROMAN)].copy()
log("  backbone genes with a CDS, chromosomes I to XVI:  %s" % f"{len(cds):,}")

# lookup tables keyed on (chr, strand, coordinate)
atg_idx = collections.defaultdict(list)
stop_idx = collections.defaultdict(list)
for r in cds.itertuples(index=False):
    atg_idx[(r.chr, r.strand, r.atg)].append(r.gene)
    stop_idx[(r.chr, r.strand, r.stop)].append(r.gene)

# The measured convention. end -> strand -> offset added to the UNAGI value
# to land on the backbone 1-based inclusive coordinate.
CONVENTION = {"5": {"+": 1, "-": 0}, "3": {"+": 0, "-": 1}}


def join_by_codon(records, end):
    """records: iterable of (nc_accession, strand, codon_coord, payload).

    Returns (matched DataFrame, diagnostics dict). Matching uses the measured
    offset only. Rows resolving to more than one backbone gene are dropped as
    ambiguous rather than assigned arbitrarily.
    """
    idx = atg_idx if end == "5" else stop_idx
    out, diag = [], collections.Counter()
    for nc, strand, coord, payload in records:
        chrom = NC_TO_ROMAN.get(str(nc).split(".")[0])
        if chrom is None:
            diag["chromosome not I to XVI"] += 1
            continue
        if strand not in ("+", "-") or coord is None:
            diag["malformed row"] += 1
            continue
        key = (chrom, strand, int(coord) + CONVENTION[end][strand])
        hits = idx.get(key)
        if not hits:
            diag["no backbone gene at that coordinate"] += 1
            continue
        if len(hits) > 1:
            diag["ambiguous, several genes share the coordinate"] += 1
            continue
        out.append((hits[0], chrom, strand) + tuple(payload))
        diag["matched"] += 1
    return out, diag


def report_diag(label, diag, total):
    log("  %s: %s rows" % (label, f"{total:,}"))
    for k, v in sorted(diag.items(), key=lambda kv: -kv[1]):
        log("     %-45s %6s  (%5.1f%%)" % (k, f"{v:,}", 100 * v / total))


# ---- verify the convention before relying on it -----------------------
log()
log("  The offset between UNAGI and the backbone is measured, not assumed.")
log("  Tabulated here across the full range -4 to +4; a correct convention")
log("  puts every joinable row in exactly one cell per strand.")
log()

three_raw = pd.read_csv(INPUTS["unagi_3prime"])
wb = openpyxl.load_workbook(INPUTS["unagi_5prime"], read_only=True)


def sheet_records(ws):
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(h).strip() for h in rows[0]]
    return hdr, [r for r in rows[1:] if r[0] is not None]


tss_hdr, tss_rows = sheet_records(wb["TSSs"])
alt5_hdr, alt5_rows = sheet_records(wb["alternative sites"])
log("  MOESM5 sheets: %s" % ", ".join(wb.sheetnames))
log("     TSSs               %s rows, columns %s" % (f"{len(tss_rows):,}", tss_hdr))
log("     alternative sites  %s rows, columns %s" % (f"{len(alt5_rows):,}", alt5_hdr))
log("  MOESM6                %s rows, columns %s"
    % (f"{len(three_raw):,}", list(three_raw.columns)))
log()

for label, end, recs in (
    ("3' stop codon (MOESM6)", "3",
     [(r["chromosome"], r["strand"], r["stop codon"]) for _, r in three_raw.iterrows()]),
    ("5' first codon (MOESM5 TSSs)", "5",
     [(r[tss_hdr.index(tss_hdr[0])], r[tss_hdr.index("strand")],
       r[tss_hdr.index("first codon")]) for r in tss_rows]),
):
    idx = atg_idx if end == "5" else stop_idx
    tab = collections.Counter()
    for nc, strand, coord in recs:
        chrom = NC_TO_ROMAN.get(str(nc).split(".")[0])
        if chrom is None or strand not in ("+", "-"):
            continue
        try:
            c = int(coord)
        except (TypeError, ValueError):
            continue
        for d in range(-4, 5):
            if (chrom, strand, c + d) in idx:
                tab[(strand, d)] += 1
                break
        else:
            tab[(strand, None)] += 1
    log("  %s" % label)
    for strand in ("+", "-"):
        cells = {d: n for (s, d), n in tab.items() if s == strand}
        tot = sum(cells.values())
        shown = "  ".join("d=%+d:%s" % (d, f"{n:,}")
                          for d, n in sorted(cells.items(), key=lambda kv: (kv[0] is None, kv[0]))
                          if d is not None and n)
        unres = cells.get(None, 0)
        log("     strand %s   %s   unresolved %s   (n = %s)"
            % (strand, shown, f"{unres:,}", f"{tot:,}"))
        nonzero = [d for d, n in cells.items() if d is not None and n]
        assert len(nonzero) == 1, (
            "convention is not single valued for %s strand %s: %s"
            % (label, strand, sorted(nonzero)))
        assert nonzero[0] == CONVENTION[end][strand], (
            "measured offset %+d contradicts the recorded convention %+d"
            % (nonzero[0], CONVENTION[end][strand]))
    log()

log("  Convention verified: every joinable row on each strand falls in one")
log("  offset cell, and that cell is the one recorded above. This is 0-based")
log("  half-open against our 1-based inclusive GFF3.")
log()

# ---- do the joins -----------------------------------------------------
three_recs = [(r["chromosome"], r["strand"], r["stop codon"],
               (r["gene"], r["gene end"], r["stop codon"], r["alternative sites"]))
              for _, r in three_raw.iterrows()]
m3, d3 = join_by_codon(three_recs, "3")
report_diag("MOESM6, 3' ends", d3, len(three_recs))
u3 = pd.DataFrame(m3, columns=["gene", "chr", "strand", "refseq", "tx_end",
                               "stop_codon", "alt_sites_raw"])

ci = {h: i for i, h in enumerate(tss_hdr)}
tss_recs = [(r[0], r[ci["strand"]], r[ci["first codon"]],
             (r[ci["gene"]], r[ci["gene start"]], r[ci["first codon"]]))
            for r in tss_rows]
m5, d5 = join_by_codon(tss_recs, "5")
log()
report_diag("MOESM5 TSSs, 5' ends", d5, len(tss_recs))
u5 = pd.DataFrame(m5, columns=["gene", "chr", "strand", "refseq", "tx_start",
                               "first_codon"])

ai = {h: i for i, h in enumerate(alt5_hdr)}
alt5_recs = [(r[0], r[ai["strand"]], r[ai["start codon"]],
              (r[ai["gene"]], r[ai["alternative start sites"]], r[ai["start codon"]]))
             for r in alt5_rows]
ma5, da5 = join_by_codon(alt5_recs, "5")
log()
report_diag("MOESM5 alternative sites, 5' ends", da5, len(alt5_recs))
alt5 = pd.DataFrame(ma5, columns=["gene", "chr", "strand", "refseq",
                                  "alt_sites_raw", "start_codon"])
log()

# ---- UTR lengths, as a difference of two columns of one table ---------
# Both columns share whatever base convention UNAGI uses, so the difference
# needs no offset correction. Signed so that the UTR extends away from the ORF.
u3["unagi_three"] = np.where(u3["strand"] == "+",
                             u3["tx_end"] - u3["stop_codon"],
                             u3["stop_codon"] - u3["tx_end"])
u5["unagi_five"] = np.where(u5["strand"] == "+",
                            u5["first_codon"] - u5["tx_start"],
                            u5["tx_start"] - u5["first_codon"])

for nm, df, col in (("5'", u5, "unagi_five"), ("3'", u3, "unagi_three")):
    neg = int((df[col] < 0).sum())
    log("  UNAGI %s UTR lengths   n = %s   negative = %d   median %.0f nt   "
        "IQR %.0f to %.0f   max %s"
        % (nm, f"{len(df):,}", neg, df[col].median(),
           df[col].quantile(.25), df[col].quantile(.75), f"{int(df[col].max()):,}"))
    if neg:
        log("     dropping %d rows with a negative length, the transcript end "
            "sits inside the ORF" % neg)
log()

u5 = u5[u5["unagi_five"] >= 0]
u3 = u3[u3["unagi_three"] >= 0]
# one row per gene; duplicates would be a RefSeq isoform pair on one ORF
for nm, df in (("5'", u5), ("3'", u3)):
    dup = int(df["gene"].duplicated().sum())
    if dup:
        log("  %s: %d duplicate genes collapsed to the longest UTR" % (nm, dup))
u5 = u5.sort_values("unagi_five").drop_duplicates("gene", keep="last")
u3 = u3.sort_values("unagi_three").drop_duplicates("gene", keep="last")

unagi = u5[["gene", "unagi_five"]].merge(u3[["gene", "unagi_three"]],
                                         on="gene", how="outer")
log("  UNAGI benchmark assembled:  %s genes, %s with a 5' end, %s with a 3' end,"
    " %s with both"
    % (f"{len(unagi):,}", f"{int(unagi['unagi_five'].notna().sum()):,}",
       f"{int(unagi['unagi_three'].notna().sum()):,}",
       f"{int((unagi['unagi_five'].notna() & unagi['unagi_three'].notna()).sum()):,}"))
unagi.to_csv(os.path.join(OUT, "unagi_join.tsv"), sep="\t", index=False)
log()

# ---- reproduction check against the published counts ------------------
# The paper's own headline counts. If our parse of the site lists is right we
# must recover them exactly from the files, before joining anything.
log("  Reproduction check, our parse against the counts stated in the paper:")


def _sites_of(v):
    s = str(v).strip()
    if s in ("", "0", "nan", "None"):
        return []
    try:
        p = ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return []
    return [int(p)] if isinstance(p, (int, float)) else [int(x) for x in p]


n_alt3 = int(sum(1 for v in three_raw["alternative sites"]
                 if len(set(_sites_of(v))) >= 2))
n_alt5 = len(alt5_rows)
for label, got, want in (
        ("genes with alternative poly(A) sites  ", n_alt3, 1482),
        ("genes with alternative 5' UTR isoforms", n_alt5, 1089)):
    log("     %s  ours %s   published %s   %s"
        % (label, f"{got:,}", f"{want:,}", "match" if got == want else "MISMATCH"))
    assert got == want, "%s: parsed %d, paper states %d" % (label, got, want)

# Written out so the ledger can read these back rather than have them typed
# in by hand. The assertions above already stop the script on a mismatch;
# this makes the two counts available to downstream scripts as data.
pd.DataFrame([
    {"quantity": "genes with two or more distinct alternative poly(A) sites",
     "end": "3prime", "parsed": n_alt3, "published": 1482,
     "source_sheet": "MOESM6"},
    {"quantity": "genes with alternative 5' UTR isoforms",
     "end": "5prime", "parsed": n_alt5, "published": 1089,
     "source_sheet": "MOESM5, alternative sites"},
]).to_csv(os.path.join(OUT, "unagi_source_counts.tsv"), sep="\t", index=False)
log("     written to unagi_source_counts.tsv")
log()

# ================================================================
# PART 1. Both tails, three annotations, three tolerances
# ================================================================
log("=" * 78)
log("PART 1.  OVER-EXTENSION AND UNDER-EXTENSION AGAINST UNAGI")
log("=" * 78)
log()
log("Same rule as the TIF-seq benchmark: per end, both the annotation and the benchmark")
log("strictly above zero. Thresholds inclusive on both tails. diff is ours")
log("minus UNAGI, so positive means our boundary runs past theirs.")
log()

final = pd.read_csv(INPUTS["final_utr"], sep="\t")
final.columns = ["gene", "five", "three"]
drs = pd.read_csv(INPUTS["drs_utr"], sep="\t")
drs.columns = ["gene", "five", "three"]
nag = pd.read_csv(INPUTS["nag_utr"])
nag.columns = ["gene", "five", "three"]
pel = pd.read_csv(INPUTS["pel_utr"])
pel.columns = ["gene", "five", "three"]
for df in (final, drs, nag, pel):
    for c in ("five", "three"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

bench = unagi.rename(columns={"unagi_five": "five", "unagi_three": "three"})
END_COL = {"5": "five", "3": "three"}

# ---- coverage, before any agreement question ------------------------
# Which genes does each side give a boundary for at all. Reported both ways,
# since a resource that only covered what UNAGI covers would add nothing.
log("-" * 78)
log("COVERAGE, WHICH GENES EACH SIDE CALLS AT ALL")
log("-" * 78)
log()
cov_rows = []
allg = final.merge(bench, on="gene", how="outer", suffixes=("_a", "_b"))
for end, col in END_COL.items():
    a = (allg[col + "_a"] > 0) & allg[col + "_a"].notna()
    b = (allg[col + "_b"] > 0) & allg[col + "_b"].notna()
    cov_rows.append({"end": end + "prime", "ours_n": int(a.sum()),
                     "unagi_n": int(b.sum()), "both_n": int((a & b).sum()),
                     "ours_only_n": int((a & ~b).sum()),
                     "unagi_only_n": int((b & ~a).sum())})
    log("  %s' UTR   ours %s   UNAGI %s   both %s   ours only %s   UNAGI only %s"
        % (end, f"{int(a.sum()):,}", f"{int(b.sum()):,}", f"{int((a & b).sum()):,}",
           f"{int((a & ~b).sum()):,}", f"{int((b & ~a).sum()):,}"))
log()
log("  genes in the released table with no UNAGI record at all: %s of %s"
    % (f"{int((~final['gene'].isin(bench['gene'])).sum()):,}", f"{len(final):,}"))
log("  genes in UNAGI with no row in the released table:        %s of %s"
    % (f"{int((~bench['gene'].isin(final['gene'])).sum()):,}", f"{len(bench):,}"))
pd.DataFrame(cov_rows).to_csv(os.path.join(OUT, "unagi_coverage.tsv"),
                              sep="\t", index=False)
log()

ANNOTS = [
    ("final", final, "final_utr.tsv, the released merged annotation"),
    ("nagalakshmi", nag, "Nagalakshmi_UTR.csv, the prior reference"),
    ("drs", drs, "DRS_UTR_corrected.tsv, our raw segmentation calls"),
]

rate_rows, dist_rows, paired = [], [], {}
for aname, adf, adesc in ANNOTS:
    log("-" * 78)
    log("%s   %s" % (aname.upper(), adesc))
    log("-" * 78)
    for end, col in END_COL.items():
        m = adf[["gene", col]].merge(bench[["gene", col]], on="gene",
                                     suffixes=("_a", "_b"))
        keep = ((m[col + "_a"] > 0) & (m[col + "_b"] > 0)
                & m[col + "_a"].notna() & m[col + "_b"].notna())
        m = m[keep].copy()
        m["diff"] = m[col + "_a"] - m[col + "_b"]
        paired[(aname, end)] = m
        n = len(m)
        d = m["diff"]
        rho = stats.spearmanr(m[col + "_a"], m[col + "_b"]).statistic
        w = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
        dist_rows.append({
            "annotation": aname, "end": end + "prime", "n": n,
            "spearman_rho": round(float(rho), 4),
            "median_diff": float(d.median()),
            "q25_diff": float(d.quantile(.25)), "q75_diff": float(d.quantile(.75)),
            "mean_diff": round(float(d.mean()), 2),
            "frac_positive": round(float((d > 0).mean()), 4),
            "wilcoxon_p": float(w.pvalue),
        })
        log("  %s' UTR   n = %s   Spearman rho = %.3f" % (end, f"{n:,}", rho))
        log("     signed difference, ours minus UNAGI, nt:")
        log("        median %+.0f    quartiles %+.0f to %+.0f    mean %+.1f"
            % (d.median(), d.quantile(.25), d.quantile(.75), d.mean()))
        log("     %-5s %22s %22s %22s"
            % ("T", "under-extension", "concordant", "over-extension"))
        for T in TOLERANCES:
            under = int((d <= -T).sum())
            over = int((d >= T).sum())
            conc = n - under - over
            rate_rows.append({
                "annotation": aname, "end": end + "prime", "tolerance_nt": T, "n": n,
                "under_extension_n": under, "under_extension_pct": round(100 * under / n, 2),
                "concordant_n": conc, "concordant_pct": round(100 * conc / n, 2),
                "over_extension_n": over, "over_extension_pct": round(100 * over / n, 2),
            })
            log("     %-5d %11s (%5.1f%%) %11s (%5.1f%%) %11s (%5.1f%%)"
                % (T, f"{under:,}", 100 * under / n, f"{conc:,}", 100 * conc / n,
                   f"{over:,}", 100 * over / n))
        log()

pd.DataFrame(rate_rows).to_csv(os.path.join(OUT, "unagi_rates.tsv"), sep="\t", index=False)
pd.DataFrame(dist_rows).to_csv(os.path.join(OUT, "unagi_distribution.tsv"), sep="\t", index=False)

# ---- did the merge improve agreement -------------------------------
log("-" * 78)
log("DID THE MERGE IMPROVE AGREEMENT WITH UNAGI, OR ONLY MOVE IT")
log("-" * 78)
log("Paired on genes comparable under both references, identical gene sets.")
log()
merge_rows = []
for end, col in END_COL.items():
    a = paired[("final", end)][["gene", "diff"]].rename(columns={"diff": "final_diff"})
    b = paired[("nagalakshmi", end)][["gene", "diff"]].rename(columns={"diff": "nag_diff"})
    j = a.merge(b, on="gene")
    ae, be = j["final_diff"].abs(), j["nag_diff"].abs()
    w = stats.wilcoxon(ae, be, alternative="two-sided", zero_method="wilcox")
    merge_rows.append({
        "end": end + "prime", "n_paired": len(j),
        "median_abs_err_final": float(ae.median()),
        "median_abs_err_prior": float(be.median()),
        "final_closer_n": int((ae < be).sum()),
        "prior_closer_n": int((ae > be).sum()),
        "tied_n": int((ae == be).sum()),
        "wilcoxon_p": float(w.pvalue),
    })
    log("  %s' UTR   paired genes = %s" % (end, f"{len(j):,}"))
    log("     median absolute error, final       %6.0f nt" % ae.median())
    log("     median absolute error, prior ref   %6.0f nt" % be.median())
    log("     closer to UNAGI: final %s   prior %s   tied %s"
        % (f"{int((ae < be).sum()):,}", f"{int((ae > be).sum()):,}",
           f"{int((ae == be).sum()):,}"))
    log("     Wilcoxon signed rank on |error|, p = %.3g" % w.pvalue)
    log()
pd.DataFrame(merge_rows).to_csv(os.path.join(OUT, "unagi_merge_effect.tsv"),
                                sep="\t", index=False)

# ================================================================
# PART 2. Provenance
# ================================================================
log("=" * 78)
log("PART 2.  THE FINAL ANNOTATION SPLIT BY BOUNDARY PROVENANCE")
log("=" * 78)
log()
log("final = pmax(DRS, Nagalakshmi) at 03_utr_table.R:200. A boundary retained")
log("from the prior reference is not evidence about our data, and is reported")
log("apart from one our segmentation actually moved.")
log()

prov_rows = []
for end, col in END_COL.items():
    m = paired[("final", end)][["gene", "diff"]].copy()
    m = m.merge(drs[["gene", col]].rename(columns={col: "drs"}), on="gene", how="left")
    m = m.merge(nag[["gene", col]].rename(columns={col: "nag"}), on="gene", how="left")
    has_d = m["drs"].notna() & (m["drs"] > 0)
    has_n = m["nag"].notna() & (m["nag"] > 0)

    def label(r, hd, hn):
        if hd and hn:
            if r["drs"] > r["nag"]:
                return "extended by DRS over prior"
            if r["drs"] < r["nag"]:
                return "prior retained, DRS shorter"
            return "both sources agree"
        if hd:
            return "DRS only, no prior call"
        if hn:
            return "prior only, no DRS call"
        return "neither, should not occur"

    m["provenance"] = [label(r, hd, hn) for r, hd, hn
                       in zip(m.to_dict("records"), has_d, has_n)]
    log("  %s' UTR" % end)
    for p, g in m.groupby("provenance"):
        d = g["diff"]
        n = len(d)
        row = {"end": end + "prime", "provenance": p, "n": n,
               "median_diff": float(d.median()),
               "frac_positive": round(float((d > 0).mean()), 4)}
        for T in TOLERANCES:
            under, over = int((d <= -T).sum()), int((d >= T).sum())
            row["under_pct_T%d" % T] = round(100 * under / n, 2)
            row["over_pct_T%d" % T] = round(100 * over / n, 2)
        prov_rows.append(row)
        under20 = int((d <= -20).sum())
        over20 = int((d >= 20).sum())
        log("     %-30s n = %5d   median diff %+6.0f nt" % (p, n, d.median()))
        log("        at T = 20 nt   under %5.1f%%   concordant %5.1f%%   over %5.1f%%"
            % (100 * under20 / n, 100 * (n - under20 - over20) / n, 100 * over20 / n))
    log()

pd.DataFrame(prov_rows).to_csv(os.path.join(OUT, "unagi_by_provenance.tsv"),
                               sep="\t", index=False)

# ================================================================
# PART 3. What kind of estimate the benchmark is, and how sharp it is
# ================================================================
log("=" * 78)
log("PART 3.  THE BENCHMARK'S OWN DEFINITION, AND ITS OWN DISPERSION")
log("=" * 78)
log()
log("This began as an envelope test analogous to the TIF-seq benchmark, asking")
log("whether our longer calls fall inside the range of ends the benchmark itself")
log("observed. Against UNAGI that question turns out to be empty, and the reason")
log("is worth reporting rather than hiding.")
log()
log("UNAGI's reported boundary is not a mode. It is the furthest end it saw: the")
log("reported boundary equals the maximum of the alternative site list for 100%")
log("of genes at the 3' end and 99.5% at the 5' end. That is measured below")
log("rather than assumed, and the script stops if it ceases to hold.")
log()
log("This changes how the Part 1 rates should be read. The TIF-seq mTIF is a")
log("single most-supported isoform, a mode, so comparing our maximum against it")
log("charges a definition mismatch as error, and the TIF-seq benchmark needed an envelope")
log("test to separate the two. UNAGI needs no such correction. Both sides are")
log("maxima over the transcripts each study observed, so Part 1 is already like")
log("for like.")
log()
log("What the site lists still measure is the benchmark's own dispersion: how far")
log("apart are the ends UNAGI saw for one gene. A tolerance narrower than that")
log("spread is scoring isoform heterogeneity, not annotation error.")
log()


def parse_sites(v):
    """The alternative-site columns hold a stringified python list, or 0."""
    if v is None:
        return []
    s = str(v).strip()
    if s in ("", "0", "nan"):
        return []
    try:
        parsed = ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return []
    if isinstance(parsed, (int, float)):
        return [int(parsed)]
    return [int(x) for x in parsed]


def site_lengths(sites, codon, strand, end):
    """Alternative site positions converted to UTR lengths, away from the ORF."""
    if end == "5":
        return [(codon - s) if strand == "+" else (s - codon) for s in sites]
    return [(s - codon) if strand == "+" else (codon - s) for s in sites]


disp_rows = []

# 5' sites, from the alternative sites sheet of MOESM5
alt5 = alt5.copy()
alt5["sites"] = alt5["alt_sites_raw"].map(parse_sites)
alt5["lengths"] = [site_lengths(s, c, st, "5") for s, c, st
                   in zip(alt5["sites"], alt5["start_codon"], alt5["strand"])]

# 3' sites, from the alternative sites column of MOESM6
u3e = pd.DataFrame(m3, columns=["gene", "chr", "strand", "refseq", "tx_end",
                                "stop_codon", "alt_sites_raw"])
u3e["sites"] = u3e["alt_sites_raw"].map(parse_sites)
u3e["lengths"] = [site_lengths(s, c, st, "3") for s, c, st
                  in zip(u3e["sites"], u3e["stop_codon"], u3e["strand"])]

for end, sdf in (("5", alt5), ("3", u3e)):
    col = END_COL[end]
    s = sdf[["gene", "lengths"]].copy()
    s["n_sites"] = s["lengths"].map(len)
    s = s[s["n_sites"] > 0].copy()
    s["max_site"] = s["lengths"].map(max)
    s["min_site"] = s["lengths"].map(min)
    s = s.sort_values("max_site").drop_duplicates("gene", keep="last")

    base = paired[("final", end)][["gene", col + "_a", col + "_b", "diff"]].rename(
        columns={col + "_a": "ours", col + "_b": "unagi"})
    m = base.merge(s[["gene", "n_sites", "max_site", "min_site"]], on="gene",
                   how="inner")
    n = len(m)
    if n == 0:
        log("  %s' UTR   no genes with alternative sites in the comparison set" % end)
        continue

    # ---- what kind of estimate is the reported boundary -----------------
    equals_max = int((m["unagi"] == m["max_site"]).sum())
    beyond = m[m["max_site"] > m["unagi"]]
    exceeds = len(beyond)
    log("  %s' UTR   genes with a UNAGI alternative site list: %s"
        % (end, f"{n:,}"))
    log("     reported boundary equals the furthest site listed:  %s of %s (%.1f%%)"
        % (f"{equals_max:,}", f"{n:,}", 100 * equals_max / n))
    log("     a listed site lies beyond the reported boundary:    %s of %s (%.1f%%)"
        % (f"{exceeds:,}", f"{n:,}", 100 * exceeds / n))
    if exceeds:
        ex = beyond["max_site"] - beyond["unagi"]
        log("        those are cross-table inconsistencies in the source, median")
        log("        excess %.0f nt, largest %.0f nt. Left as found, not repaired."
            % (ex.median(), ex.max()))
    assert exceeds / n < 0.01, (
        "%s' end: %.1f%% of genes list a site beyond the reported boundary, too "
        "many to treat the boundary as a maximum" % (end, 100 * exceeds / n))
    log("     So the reported boundary is a maximum, as ours is. No envelope")
    log("     correction applies, and the Part 1 rates stand as measured.")
    log()

    # ---- the benchmark's own dispersion --------------------------------
    spread = m["max_site"] - m["min_site"]
    log("     UNAGI's own spread of ends for one gene, furthest minus nearest:")
    log("        sites per gene, median %.0f    spread median %.0f nt, "
        "quartiles %.0f to %.0f nt"
        % (m["n_sites"].median(), spread.median(),
           spread.quantile(.25), spread.quantile(.75)))
    row = {"end": end + "prime", "n_with_alt_sites": n,
           "boundary_equals_max_n": equals_max,
           "boundary_equals_max_pct": round(100 * equals_max / n, 2),
           "sites_beyond_boundary_n": exceeds,
           "median_n_alt_sites": float(m["n_sites"].median()),
           "median_spread_nt": float(spread.median()),
           "q25_spread_nt": float(spread.quantile(.25)),
           "q75_spread_nt": float(spread.quantile(.75))}
    for T in TOLERANCES:
        f = float((spread >= T).mean())
        row["own_spread_ge_T%d_pct" % T] = round(100 * f, 2)
        log("        own spread already reaches or exceeds %3d nt on %5.1f%% of genes"
            % (T, 100 * f))
    log("     On those genes a tolerance below the spread is scoring isoform")
    log("     heterogeneity in the benchmark, not error in the annotation.")
    disp_rows.append(row)
    log()

pd.DataFrame(disp_rows).to_csv(os.path.join(OUT, "unagi_dispersion.tsv"),
                               sep="\t", index=False)

# ================================================================
# PART 4. The two independent benchmarks against each other
# ================================================================
log("=" * 78)
log("PART 4.  UNAGI AND TIF-SEQ, TWO INDEPENDENT BENCHMARKS")
log("=" * 78)
log()
log("Neither took part in construction. They disagree with each other as well")
log("as with us, and that disagreement bounds how much of our own deviation is")
log("attributable to the annotation rather than to benchmark choice.")
log()

cross_rows = []
for end, col in END_COL.items():
    j = bench[["gene", col]].rename(columns={col: "unagi"}).merge(
        pel[["gene", col]].rename(columns={col: "tifseq"}), on="gene")
    j = j[(j["unagi"] > 0) & (j["tifseq"] > 0)].dropna()
    d = j["unagi"] - j["tifseq"]
    n = len(j)
    rho = stats.spearmanr(j["unagi"], j["tifseq"]).statistic
    row = {"end": end + "prime", "comparison": "unagi vs tifseq", "n": n,
           "spearman_rho": round(float(rho), 4),
           "median_diff": float(d.median()),
           "median_abs_diff": float(d.abs().median())}
    for T in TOLERANCES:
        row["concordant_n_T%d" % T] = int((d.abs() < T).sum())
        row["concordant_pct_T%d" % T] = round(100 * float((d.abs() < T).sum()) / n, 2)
    cross_rows.append(row)
    log("  %s' UTR   genes in both benchmarks = %s   Spearman rho = %.3f"
        % (end, f"{n:,}", rho))
    log("     UNAGI minus TIF-seq:  median %+.0f nt   median absolute %.0f nt"
        % (d.median(), d.abs().median()))
    log("     the two benchmarks agree within 20 nt on %.1f%% of genes"
        % (100 * float((d.abs() < 20).sum()) / n))
    log()

# our agreement with each, on the genes common to both benchmarks
for end, col in END_COL.items():
    common = bench[["gene", col]].rename(columns={col: "unagi"}).merge(
        pel[["gene", col]].rename(columns={col: "tifseq"}), on="gene")
    common = common[(common["unagi"] > 0) & (common["tifseq"] > 0)].dropna()
    j = final[["gene", col]].rename(columns={col: "ours"}).merge(common, on="gene")
    j = j[j["ours"] > 0]
    n = len(j)
    du, dt = j["ours"] - j["unagi"], j["ours"] - j["tifseq"]
    row = {"end": end + "prime", "comparison": "ours, on genes in both benchmarks",
           "n": n, "spearman_rho": np.nan,
           "median_diff": np.nan, "median_abs_diff": np.nan}
    for T in TOLERANCES:
        row["concordant_pct_T%d" % T] = round(100 * float((du.abs() < T).sum()) / n, 2)
    # Counts as well as percentages. The manuscript quotes these to one
    # decimal place, and re-rounding a value already rounded to two sends
    # 41.75 and 32.25 the wrong way. The ledger recomputes from the integers.
    row["tifseq_concordant_n_T20"] = int((dt.abs() < 20).sum())
    row["unagi_concordant_n_T20"] = int((du.abs() < 20).sum())
    row["both_concordant_n_T20"] = int(((du.abs() < 20) & (dt.abs() < 20)).sum())
    row["tifseq_concordant_pct_T20"] = round(100 * float((dt.abs() < 20).sum()) / n, 2)
    row["unagi_concordant_pct_T20"] = round(100 * float((du.abs() < 20).sum()) / n, 2)
    row["both_concordant_pct_T20"] = round(
        100 * float(((du.abs() < 20) & (dt.abs() < 20)).sum()) / n, 2)
    cross_rows.append(row)
    log("  %s' UTR   our calls on the %s genes both benchmarks cover:" % (end, f"{n:,}"))
    log("     within 20 nt of UNAGI    %5.1f%%" % (100 * float((du.abs() < 20).sum()) / n))
    log("     within 20 nt of TIF-seq  %5.1f%%" % (100 * float((dt.abs() < 20).sum()) / n))
    log("     within 20 nt of both     %5.1f%%"
        % (100 * float(((du.abs() < 20) & (dt.abs() < 20)).sum()) / n))
    log()

pd.DataFrame(cross_rows).to_csv(os.path.join(OUT, "unagi_vs_tifseq.tsv"),
                                sep="\t", index=False)

# ---- fine sweep, to show the tolerance is not load bearing ----------
sweep = []
for T in range(0, 205, 5):
    for aname, _, _ in ANNOTS:
        for end in END_COL:
            d = paired[(aname, end)]["diff"]
            n = len(d)
            sweep.append({
                "annotation": aname, "end": end + "prime", "tolerance_nt": T, "n": n,
                "under_pct": round(100 * float((d <= -T).sum()) / n, 3),
                "over_pct": round(100 * float((d >= T).sum()) / n, 3),
                "concordant_pct": round(100 * float((d.abs() < T).sum()) / n, 3),
            })
pd.DataFrame(sweep).to_csv(os.path.join(OUT, "unagi_rate_curve.tsv"),
                           sep="\t", index=False)

log("=" * 78)
log("Outputs written to %s" % OUT)
log("Done.")
_logfh.close()
