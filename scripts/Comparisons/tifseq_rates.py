#!/usr/bin/env python3
"""
tifseq_rates.py
==================
The TIF-seq validation reported with explicit over-extension AND
under-extension rates (not just one tail), for three annotations, at
three tolerances, split by the provenance of each boundary and by the
read support of the benchmark itself. Underlies Supplementary Fig. S6
(plotted by plot_suppl_fig6.py).

WHAT IS COMPARED
------------------
  final    final_utr.tsv            the released merged annotation
  nagalakshmi  Nagalakshmi_UTR.csv  the prior reference, one side of the merge
  drs      DRS_UTR_corrected.tsv    our raw segmentation calls, the other side
  benchmark  Pelechano_UTR.csv (+ .qc.csv sidecar)  TIF-seq major covering
             isoform per ORF, derived by derive_pelechano_utr.py, used in
             construction by nothing

Reporting the prior reference against the same benchmark is the
comparison that answers whether the merge improved agreement with
TIF-seq or only moved it.

THE COMPARISON SET
---------------------
Per end, both the annotation and TIF-seq strictly greater than zero
(not "both sides"), matching plot_final_vs_tifseq_comparison.R.
Categories at tolerance T, on diff = ours minus TIF-seq:
  diff <= -T   under-extension, TIF-seq longer, we truncate
  |diff| < T   concordant
  diff >= +T   over-extension, TIF-seq shorter, we run past it
The threshold is inclusive on both tails.

THE TOLERANCE, DERIVED NOT ASSERTED
--------------------------------------
TIF-seq resolution is not a protocol constant to cite, so it's measured
from the raw TIF-seq table: for each gene, every isoform spanning the
ORF is weighted by its read count, and we ask how far those reads place
the end from the end of that gene's own major isoform (Part 1). A
tolerance narrower than that dispersion measures TIF-seq isoform
heterogeneity, not annotation error.

Input:  Data/Annotation/{final_utr.tsv, DRS_UTR_corrected.tsv,
        Nagalakshmi_UTR.csv, Pelechano_UTR.csv, Pelechano_UTR.qc.csv,
        gene_models_backbone_with_UTR_slots.gff3}
        the raw TIF-seq isoform table GSE39128_tsedall.txt.gz (GEO
        GSE39128), the same file derive_pelechano_utr.py consumes
Output: tifseq_rates.tsv            main result, annotation x end x tolerance
        tifseq_distribution.tsv     signed difference distributions
        tifseq_by_provenance.tsv    final annotation split by boundary source
        tifseq_by_support.tsv       stratified by mTIF read support
        tifseq_resolution.tsv       measured dispersion of the benchmark
        tifseq_resolution_curve.tsv resolution as a function of tolerance
        tifseq_rate_curve.tsv       fine tolerance sweep, for plotting
        tifseq_envelope.tsv         envelope test (Part 5)

Run:
  python3 tifseq_rates.py --annotation-dir Data/Annotation \\
      --pelechano-raw GSE39128_tsedall.txt.gz --outdir Data/TIFseq
"""

import argparse
import gzip
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

TOLERANCES = [20, 50, 100]
RES_GRID = list(range(0, 205, 5))
ROMAN_TO_ARABIC = {r: i + 1 for i, r in enumerate(
    ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
     "XI", "XII", "XIII", "XIV", "XV", "XVI"])}
COUNT_COLS = ["ypd", "gal", "lypd", "lgal", "nypd", "ngal"]
END_COL = {"5": "five", "3": "three"}


def log(msg=""):
    print(msg)


def load_orfs(path):
    rows = []
    id_re = re.compile(r"ID=([^;]+)")
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = ROMAN_TO_ARABIC.get(f[0].strip())
            if chrom is None:
                continue
            m = id_re.search(f[8])
            if not m:
                continue
            rows.append((m.group(1), chrom, f[6], int(f[3]), int(f[4])))
    return pd.DataFrame(rows, columns=["gene", "chr", "strand", "start", "end"])


def wquantile(vals, wts, q):
    """Weighted quantile. vals and wts are 1D arrays, wts non-negative."""
    o = np.argsort(vals)
    v, w = vals[o], wts[o]
    cw = np.cumsum(w)
    if cw[-1] <= 0:
        return np.nan
    return float(np.interp(q * cw[-1], cw - 0.5 * w, v))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--pelechano-raw", required=True,
                     help="GSE39128_tsedall.txt.gz, the raw TIF-seq isoform table")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    ANN, OUT = args.annotation_dir, args.outdir
    os.makedirs(OUT, exist_ok=True)

    final = pd.read_csv(os.path.join(ANN, "final_utr.tsv"), sep="\t")
    final.columns = ["gene", "five", "three"]
    drs = pd.read_csv(os.path.join(ANN, "DRS_UTR_corrected.tsv"), sep="\t")
    drs.columns = ["gene", "five", "three"]
    nag = pd.read_csv(os.path.join(ANN, "Nagalakshmi_UTR.csv"))
    nag.columns = ["gene", "five", "three"]
    pel = pd.read_csv(os.path.join(ANN, "Pelechano_UTR.csv"))
    pel.columns = ["gene", "five", "three"]
    pel_qc = pd.read_csv(os.path.join(ANN, "Pelechano_UTR.qc.csv"))
    for df in (final, drs, nag, pel):
        for c in ("five", "three"):
            df[c] = pd.to_numeric(df[c], errors="coerce")

    log("TABLES")
    log("  final_utr.tsv          %6d rows" % len(final))
    log("  DRS_UTR_corrected.tsv  %6d rows" % len(drs))
    log("  Nagalakshmi_UTR.csv    %6d rows" % len(nag))
    log("  Pelechano_UTR.csv      %6d rows, %6d with a covering mTIF"
        % (len(pel), int(pel["five"].notna().sum())))
    log()

    # ============================================================
    # PART 1. Resolution of the benchmark, measured from the raw table
    # ============================================================
    log("PART 1.  HOW SHARP IS THE TIF-SEQ BOUNDARY ITSELF")
    orfs = load_orfs(os.path.join(ANN, "gene_models_backbone_with_UTR_slots.gff3"))
    log("  ORFs from the backbone, chromosomes 1 to 16:  %d" % len(orfs))

    with gzip.open(args.pelechano_raw, "rt") as fh:
        tifs = pd.read_csv(fh, sep="\t")
    tifs["total"] = tifs[COUNT_COLS].sum(axis=1)
    tifs = tifs[["chr", "strand", "t5", "t3", "total"]]
    log("  raw TIF-seq isoforms:  %s rows, %s reads total"
        % (f"{len(tifs):,}", f"{int(tifs['total'].sum()):,}"))

    groups = {}
    for (c, s), g in tifs.groupby(["chr", "strand"]):
        groups[(c, s)] = (g["t5"].to_numpy(), g["t3"].to_numpy(), g["total"].to_numpy())

    res_rows = []
    for row in orfs.itertuples(index=False):
        key = (row.chr, row.strand)
        if key not in groups:
            continue
        t5, t3, tot = groups[key]
        if row.strand == "+":
            mask = (t5 <= row.start) & (t3 >= row.end)
        else:
            mask = (t5 >= row.end) & (t3 <= row.start)
        if not mask.any():
            continue
        idx = np.flatnonzero(mask)
        w = tot[idx].astype(float)
        if w.sum() <= 0:
            continue
        best = idx[np.argmax(tot[idx])]
        a5, a3 = t5[idx].astype(float), t3[idx].astype(float)
        d5 = np.abs(a5 - float(t5[best]))
        d3 = np.abs(a3 - float(t3[best]))
        if row.strand == "+":
            u5, u3 = row.start - a5, a3 - row.end
        else:
            u5, u3 = a5 - row.end, row.start - a3
        r = {
            "gene": row.gene, "n_isoforms": len(idx), "n_reads": float(w.sum()),
            "mtif_support": int(tot[best]),
            "iqr5": wquantile(a5, w, 0.75) - wquantile(a5, w, 0.25),
            "iqr3": wquantile(a3, w, 0.75) - wquantile(a3, w, 0.25),
            "utr5_max": float(u5.max()), "utr3_max": float(u3.max()),
            "utr5_p95": wquantile(u5, w, 0.95), "utr3_p95": wquantile(u3, w, 0.95),
            "utr5_p99": wquantile(u5, w, 0.99), "utr3_p99": wquantile(u3, w, 0.99),
        }
        for T in RES_GRID:
            r["f5_within_%d" % T] = float(w[d5 <= T].sum() / w.sum())
            r["f3_within_%d" % T] = float(w[d3 <= T].sum() / w.sum())
        res_rows.append(r)

    resd = pd.DataFrame(res_rows)
    _keep = [c for c in resd.columns if not c.startswith(("f5_within", "f3_within"))] + \
            ["f%s_within_%d" % (e, T) for T in TOLERANCES for e in ("5", "3")]
    resd[_keep].to_csv(os.path.join(OUT, "tifseq_resolution.tsv"), sep="\t", index=False)
    log("  genes with at least one spanning isoform:  %d" % len(resd))
    log("  read-weighted interquartile range of ends within a gene:")
    log("     5' end   median %.0f nt   3' end   median %.0f nt"
        % (resd["iqr5"].median(), resd["iqr3"].median()))
    log()

    # ============================================================
    # PART 2. Both tails, three annotations, three tolerances
    # ============================================================
    log("PART 2.  OVER-EXTENSION AND UNDER-EXTENSION, BOTH REPORTED")
    ANNOTS = [("final", final), ("nagalakshmi", nag), ("drs", drs)]
    rate_rows, dist_rows, paired = [], [], {}

    for aname, adf in ANNOTS:
        for end, col in END_COL.items():
            m = adf[["gene", col]].merge(pel[["gene", col]], on="gene", suffixes=("_a", "_p"))
            keep = (m[col + "_a"] > 0) & (m[col + "_p"] > 0) & m[col + "_a"].notna() & m[col + "_p"].notna()
            m = m[keep].copy()
            m["diff"] = m[col + "_a"] - m[col + "_p"]
            paired[(aname, end)] = m
            n = len(m)
            rho = stats.spearmanr(m[col + "_a"], m[col + "_p"]).statistic
            d = m["diff"]
            w = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
            dist_rows.append({
                "annotation": aname, "end": end + "prime", "n": n,
                "spearman_rho": round(float(rho), 4), "median_diff": float(d.median()),
                "q25_diff": float(d.quantile(.25)), "q75_diff": float(d.quantile(.75)),
                "mean_diff": round(float(d.mean()), 2),
                "frac_positive": round(float((d > 0).mean()), 4), "wilcoxon_p": float(w.pvalue),
            })
            log("  %s %s' UTR   n = %s   Spearman rho = %.3f" % (aname, end, f"{n:,}", rho))
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
    pd.DataFrame(rate_rows).to_csv(os.path.join(OUT, "tifseq_rates.tsv"), sep="\t", index=False)
    pd.DataFrame(dist_rows).to_csv(os.path.join(OUT, "tifseq_distribution.tsv"), sep="\t", index=False)
    log()

    log("DID THE MERGE IMPROVE AGREEMENT, OR ONLY MOVE IT")
    for end, col in END_COL.items():
        a = paired[("final", end)][["gene", col + "_a", "diff"]].rename(
            columns={col + "_a": "final_len", "diff": "final_diff"})
        b = paired[("nagalakshmi", end)][["gene", col + "_a", "diff"]].rename(
            columns={col + "_a": "nag_len", "diff": "nag_diff"})
        j = a.merge(b, on="gene")
        ae, be = j["final_diff"].abs(), j["nag_diff"].abs()
        w = stats.wilcoxon(ae, be, alternative="two-sided", zero_method="wilcox")
        log("  %s' UTR  paired=%d  median|err| final=%.0f prior=%.0f  closer: final=%d prior=%d tied=%d  p=%.3g"
            % (end, len(j), ae.median(), be.median(), int((ae < be).sum()),
               int((ae > be).sum()), int((ae == be).sum()), w.pvalue))
    log()

    # ============================================================
    # PART 3. Where did each boundary come from
    # ============================================================
    log("PART 3.  THE FINAL ANNOTATION SPLIT BY BOUNDARY PROVENANCE")
    prov_rows = []
    for end, col in END_COL.items():
        m = paired[("final", end)][["gene", col + "_a", col + "_p", "diff"]].copy()
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

        m["provenance"] = [label(r, hd, hn) for r, hd, hn in
                           zip(m.to_dict("records"), has_d, has_n)]
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
            log("  %s' %-28s n=%5d  median %+6.0f nt" % (end, p, n, d.median()))
    pd.DataFrame(prov_rows).to_csv(os.path.join(OUT, "tifseq_by_provenance.tsv"), sep="\t", index=False)
    log()

    # ============================================================
    # PART 4. Stratified by mTIF read support
    # ============================================================
    log("PART 4.  STRATIFIED BY MTIF READ SUPPORT")
    BINS = [(1, 1), (2, 4), (5, 19), (20, 10 ** 9)]
    sup_rows = []
    for end, col in END_COL.items():
        m = paired[("final", end)][["gene", "diff"]].merge(
            pel_qc[["Gene", "mtif_support"]].rename(columns={"Gene": "gene"}), on="gene")
        for lo, hi in BINS:
            g = m[(m["mtif_support"] >= lo) & (m["mtif_support"] <= hi)]
            n = len(g)
            if n == 0:
                continue
            d = g["diff"]
            under, over = int((d <= -20).sum()), int((d >= 20).sum())
            lab = "%d" % lo if lo == hi else ("%d+" % lo if hi > 10 ** 8 else "%d to %d" % (lo, hi))
            sup_rows.append({"end": end + "prime", "support_bin": lab, "n": n,
                             "median_diff": float(d.median()),
                             "under_pct_T20": round(100 * under / n, 2),
                             "over_pct_T20": round(100 * over / n, 2)})
    pd.DataFrame(sup_rows).to_csv(os.path.join(OUT, "tifseq_by_support.tsv"), sep="\t", index=False)
    log()

    # ============================================================
    # PART 5. The envelope test
    # ============================================================
    log("PART 5.  THE ENVELOPE TEST")
    env_rows = []
    for end, num in (("5", "utr5"), ("3", "utr3")):
        col = END_COL[end]
        m = paired[("final", end)][["gene", col + "_a", col + "_p", "diff"]].rename(
            columns={col + "_a": "ours", col + "_p": "mtif"})
        m = m.merge(resd[["gene", num + "_max", num + "_p95", num + "_p99"]], on="gene", how="inner")
        n = len(m)
        over = m[m["diff"] >= 20]
        for lab, ref in (("most extreme observed", num + "_max"),
                         ("99th percentile of reads", num + "_p99"),
                         ("95th percentile of reads", num + "_p95")):
            inside_all = int((m["ours"] <= m[ref]).sum())
            inside_over = int((over["ours"] <= over[ref]).sum())
            env_rows.append({
                "end": end + "prime", "envelope": lab, "n": n,
                "inside_n": inside_all, "inside_pct": round(100 * inside_all / n, 2),
                "n_over_extended": len(over), "over_inside_n": inside_over,
                "over_inside_pct": round(100 * inside_over / len(over), 2) if len(over) else float("nan"),
            })
        log("  %s' UTR  n=%d  beyond every observed isoform: %d (%.1f%%)"
            % (end, n, len(m[m["ours"] > m[num + "_max"]]),
               100 * len(m[m["ours"] > m[num + "_max"]]) / n))
    pd.DataFrame(env_rows).to_csv(os.path.join(OUT, "tifseq_envelope.tsv"), sep="\t", index=False)

    # ---- fine sweeps, for plotting ----------------------------------------
    sweep = []
    for T in range(0, 205, 5):
        for aname, _ in ANNOTS:
            for end in END_COL:
                d = paired[(aname, end)]["diff"]
                n = len(d)
                sweep.append({
                    "annotation": aname, "end": end + "prime", "tolerance_nt": T, "n": n,
                    "under_pct": round(100 * float((d <= -T).sum()) / n, 3),
                    "over_pct": round(100 * float((d >= T).sum()) / n, 3),
                    "concordant_pct": round(100 * float((d.abs() < T).sum()) / n, 3),
                })
    pd.DataFrame(sweep).to_csv(os.path.join(OUT, "tifseq_rate_curve.tsv"), sep="\t", index=False)

    pd.DataFrame([{
        "tolerance_nt": T,
        "median_frac_reads_within_5prime": round(float(resd["f5_within_%d" % T].median()), 4),
        "median_frac_reads_within_3prime": round(float(resd["f3_within_%d" % T].median()), 4),
    } for T in RES_GRID]).to_csv(os.path.join(OUT, "tifseq_resolution_curve.tsv"), sep="\t", index=False)

    log()
    log("Outputs written to %s" % OUT)


if __name__ == "__main__":
    main()
