#!/usr/bin/env python3
"""
neighbour_overlap.py
======================
Quantifies how many gene pairs newly overlap as a consequence of
extending UTR boundaries, and which of those overlaps land inside a
neighbour's coding sequence -- the configuration that actually risks
confounding gene-level quantification.

WHAT COUNTS AS AN OVERLAP
--------------------------
The released reference and the prior reference are two transcriptome
FASTAs built by gffread from two GFF3 files that differ only in the
coordinates of the UTR slots. So the comparison is between those two
files directly, not between the bare ORF set and the released file:

  old span   Nagalakshmi_UTRs.gff3   ORF plus prior reference UTRs
  new span   Final_UTRs.gff3         ORF plus released UTRs

Both are read as the union of every feature belonging to a gene's mRNA
parent, which is exactly what gffread extracted. A gene pair is NEWLY
OVERLAPPING when the two spans overlap under the released annotation
and did not overlap under the prior one. Pairs that already overlapped
are reported separately, as the baseline.

Two things are counted:
  A. transcript against transcript, both members sequences in the
     released FASTA -- the multi-mapping concern in its exact form.
  B. transcript against a neighbour's ORF body, same strand -- worse
     than a UTR-into-UTR overlap, because the ORF body is where the
     reads that matter come from. This is the "worst tier" that
     abundance_concordance.py tests for a measurable quantification
     effect.

Strand is kept separate throughout: direct RNA sequencing is strand
specific, so an antisense overlap does not create the assignment
ambiguity a same-strand overlap does.

This script does not classify genes by SGD ORF status (Verified,
Uncharacterized, Dubious) or check overlap against non-mRNA features
(tRNA, snoRNA, etc.) -- both are a separate, additional analysis with
its own external SGD annotation input.

Input:  Data/Annotation/Final_UTRs.gff3, Data/Annotation/Nagalakshmi_UTRs.gff3
Output: neighbour_overlap_summary.tsv   headline counts
        neighbour_overlap_pairs.tsv     every newly overlapping pair
        neighbour_overlap_by_gene.tsv   per gene, worst overlap it gained

Run:
  python3 neighbour_overlap.py --annotation-dir Data/Annotation --outdir Data/Overlap
"""

import argparse
import os

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--annotation-dir", required=True,
                     help="Data/Annotation, holding Final_UTRs.gff3 and Nagalakshmi_UTRs.gff3")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
ANN, OUT = args.annotation_dir, args.outdir
os.makedirs(OUT, exist_ok=True)

INPUTS = {
    "final_gff": os.path.join(ANN, "Final_UTRs.gff3"),
    "nag_gff": os.path.join(ANN, "Nagalakshmi_UTRs.gff3"),
}


def log(msg=""):
    print(msg)


# ------------------------------------------------------- 1. read spans
def read_spans(path):
    """
    Per-gene transcript span, as the union of every feature belonging to
    that gene's mRNA parent. Rows with NA coordinates are absent slots
    and are skipped, matching what gffread does. The ORF body is taken
    from features whose ID ends in _CDS or _intron, not from the
    type=gene row, since both GFF3 files also label the UTR slots as
    type CDS so gffread will extract them.
    """
    lo, hi, olo, ohi, chrom, strand = {}, {}, {}, {}, {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if len(f) < 9:
                continue
            attr = f[8]
            if f[2] == "gene":
                continue
            if "Parent=" not in attr:
                continue
            par = attr.split("Parent=")[1].split(";")[0]
            if not par.endswith("_mRNA"):
                continue
            gid = par[:-5]
            chrom.setdefault(gid, f[0])
            strand.setdefault(gid, f[6])
            if f[3] == "NA" or f[4] == "NA":
                continue
            s, e = int(f[3]), int(f[4])
            lo[gid] = min(lo.get(gid, s), s)
            hi[gid] = max(hi.get(gid, e), e)
            fid = attr.split("ID=")[1].split(";")[0]
            if fid.endswith("_CDS") or fid.endswith("_intron"):
                olo[gid] = min(olo.get(gid, s), s)
                ohi[gid] = max(ohi.get(gid, e), e)
    genes = sorted(lo)
    return pd.DataFrame([{
        "Gene": g, "chrom": chrom.get(g), "strand": strand.get(g),
        "start": lo[g], "end": hi[g],
        "orf_start": olo.get(g, lo[g]), "orf_end": ohi.get(g, hi[g]),
    } for g in genes])


new = read_spans(INPUTS["final_gff"])
old = read_spans(INPUTS["nag_gff"])
log("released reference transcripts: %d" % len(new))
log("prior reference transcripts:    %d" % len(old))

spans = old[["Gene", "chrom", "strand", "orf_start", "orf_end", "start", "end"]] \
    .rename(columns={"start": "old_start", "end": "old_end"}) \
    .merge(new[["Gene", "start", "end"]].rename(columns={"start": "new_start", "end": "new_end"}),
           on="Gene", how="outer")
for c in ("old_start", "old_end", "new_start", "new_end", "orf_start", "orf_end"):
    spans[c] = spans[c].astype(int)
spans["new_len"] = spans["new_end"] - spans["new_start"] + 1
spans["old_len"] = spans["old_end"] - spans["old_start"] + 1
sp0 = spans.set_index("Gene")


# ------------------------------------------------------- 2. overlap engine
def ovl(a1, a2, b1, b2):
    return min(a2, b2) - max(a1, b1) + 1


def pairs_within(df, s_col, e_col):
    """All overlapping pairs among the spans in df, by chromosome, sweep line."""
    res = {}
    for chrom, sub in df.groupby("chrom"):
        sub = sub.sort_values(s_col)
        starts, ends, names = sub[s_col].values, sub[e_col].values, sub["Gene"].values
        active = []
        for i in range(len(sub)):
            si = starts[i]
            active = [j for j in active if ends[j] >= si]
            for j in active:
                o = ovl(starts[j], ends[j], si, ends[i])
                if o > 0:
                    res[(names[j], names[i]) if names[j] < names[i]
                        else (names[i], names[j])] = o
            active.append(i)
    return res


orf_pairs = pairs_within(spans, "orf_start", "orf_end")
old_pairs = pairs_within(spans, "old_start", "old_end")
new_pairs = pairs_within(spans, "new_start", "new_end")
newly = {k: v for k, v in new_pairs.items() if k not in old_pairs}

log("")
log("TRANSCRIPT AGAINST TRANSCRIPT")
log("  overlapping pairs, ORF bodies only      %d" % len(orf_pairs))
log("  overlapping pairs, prior reference      %d" % len(old_pairs))
log("  overlapping pairs, released reference   %d" % len(new_pairs))
log("  newly overlapping pairs                 %d" % len(newly))


def worst_tier(pair_dict, s_col, e_col):
    """Same-strand pairs whose span reaches into the partner's ORF body."""
    n = 0
    for (a, b) in pair_dict:
        ra, rb = sp0.loc[a], sp0.loc[b]
        if ra["strand"] != rb["strand"]:
            continue
        into = max(ovl(ra[s_col], ra[e_col], rb["orf_start"], rb["orf_end"]),
                   ovl(rb[s_col], rb[e_col], ra["orf_start"], ra["orf_end"]))
        if into > 0:
            n += 1
    return n


log("  same strand, into partner ORF body:")
log("    ORF bodies only      %d" % worst_tier(orf_pairs, "orf_start", "orf_end"))
log("    prior reference      %d" % worst_tier(old_pairs, "old_start", "old_end"))
log("    released reference   %d" % worst_tier(new_pairs, "new_start", "new_end"))

# ---------------------------------------------------- 3. annotate the pairs
rows = []
for (a, b), o in sorted(newly.items(), key=lambda kv: -kv[1]):
    ra, rb = sp0.loc[a], sp0.loc[b]
    same = ra["strand"] == rb["strand"]
    into_orf = max(
        ovl(ra["new_start"], ra["new_end"], rb["orf_start"], rb["orf_end"]),
        ovl(rb["new_start"], rb["new_end"], ra["orf_start"], ra["orf_end"]))
    if same:
        geom = "tandem"
    else:
        left, right = (ra, rb) if ra["orf_start"] < rb["orf_start"] else (rb, ra)
        geom = "convergent" if left["strand"] == "+" else "divergent"
    rows.append({
        "gene_a": a, "gene_b": b, "strand_a": ra["strand"], "strand_b": rb["strand"],
        "same_strand": same, "geometry": geom, "overlap_nt": int(o),
        "overlap_into_orf_nt": int(max(into_orf, 0)),
    })
pairs = pd.DataFrame(rows)
pairs.to_csv(os.path.join(OUT, "neighbour_overlap_pairs.tsv"), sep="\t", index=False)

worst = pairs[(pairs["overlap_into_orf_nt"] > 0) & pairs["same_strand"]]
log("")
log("THE NEWLY OVERLAPPING PAIRS")
log("  same strand pairs          %d" % int(pairs["same_strand"].sum()))
log("  opposite strand pairs      %d" % int((~pairs["same_strand"]).sum()))
log("  reaching into an ORF body  %d" % int((pairs["overlap_into_orf_nt"] > 0).sum()))
log("  same strand AND into ORF   %d (the worst tier)" % len(worst))

# ------------------------------------------------------- 4. per-gene view
gene_rows = []
from collections import defaultdict
inv = defaultdict(list)
for _, r in pairs.iterrows():
    inv[r["gene_a"]].append(r)
    inv[r["gene_b"]].append(r)
for g in spans["Gene"]:
    ps = inv.get(g, [])
    gene_rows.append({
        "Gene": g,
        "span_growth_nt": int(sp0.loc[g, "new_len"] - sp0.loc[g, "old_len"]),
        "n_new_gene_overlaps": len(ps),
        "max_gene_overlap_nt": int(max([p["overlap_nt"] for p in ps], default=0)),
        "max_same_strand_orf_overlap_nt": int(max(
            [p["overlap_into_orf_nt"] for p in ps if p["same_strand"]], default=0)),
    })
gdf = pd.DataFrame(gene_rows)
gdf.to_csv(os.path.join(OUT, "neighbour_overlap_by_gene.tsv"), sep="\t", index=False)

n_any = int((gdf["n_new_gene_overlaps"] > 0).sum())
log("")
log("PER GENE")
log("  genes in the released reference               %d" % len(gdf))
log("  genes that gained no new overlap               %d  (%.1f%%)"
    % (len(gdf) - n_any, 100.0 * (len(gdf) - n_any) / len(gdf)))
log("  genes that gained at least one                 %d  (%.1f%%)"
    % (n_any, 100.0 * n_any / len(gdf)))
log("  genes with a same strand overlap into an ORF   %d  (%.2f%%)"
    % (int((gdf["max_same_strand_orf_overlap_nt"] > 0).sum()),
       100.0 * (gdf["max_same_strand_orf_overlap_nt"] > 0).mean()))

summary = pd.DataFrame([
    {"measure": "gene pairs overlapping, ORF bodies only", "value": len(orf_pairs)},
    {"measure": "gene pairs overlapping, prior reference", "value": len(old_pairs)},
    {"measure": "gene pairs overlapping, released reference", "value": len(new_pairs)},
    {"measure": "newly overlapping gene pairs", "value": len(newly)},
    {"measure": "newly overlapping, same strand", "value": int(pairs["same_strand"].sum())},
    {"measure": "newly overlapping reaching an ORF body", "value": int((pairs["overlap_into_orf_nt"] > 0).sum())},
    {"measure": "newly overlapping, same strand into ORF body (worst tier)", "value": len(worst)},
    {"measure": "median new overlap length nt", "value": float(pairs["overlap_nt"].median())},
    {"measure": "genes gaining at least one new overlap", "value": n_any},
    {"measure": "genes gaining no new overlap", "value": len(gdf) - n_any},
    {"measure": "genes in released reference", "value": len(gdf)},
])
summary.to_csv(os.path.join(OUT, "neighbour_overlap_summary.tsv"), sep="\t", index=False)

log("")
log("wrote %s" % os.path.join(OUT, "neighbour_overlap_summary.tsv"))
log("wrote %s (%d pairs)" % (os.path.join(OUT, "neighbour_overlap_pairs.tsv"), len(pairs)))
log("wrote %s (%d genes)" % (os.path.join(OUT, "neighbour_overlap_by_gene.tsv"), len(gdf)))
