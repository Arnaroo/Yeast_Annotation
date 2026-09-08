#!/usr/bin/env python3
"""
read_geometry.py
===================
Read-level confirmation of the mechanism inferred in shortfall.py, with
no annotation statistic in the middle. For any primary sense alignment
in the our_ref BAM:

    start_off = POS - 1                bases between the annotated 5' end
                                       and the first aligned base
    end_off   = LN - (POS + reflen - 1) bases between the last aligned base
                                       and the annotated 3' end

Direct RNA sequencing is read 3' to 5', so end_off measures where
sequencing started and start_off measures how far it got. Processivity
loss, if it exists, is visible as start_off inflating on long
transcripts while end_off does not -- measured directly from the
alignments, no annotation quantity involved.

Not circular with the calls: this BAM should be the independent
validation library (e.g. SRR32918129), not one of the construction
libraries that produced the segmentation calls.

Input:  a BAM aligned to the our_ref transcriptome (validation library)
        shortfall_genes.tsv (shortfall.py)
        samtools on PATH
Output: read_geometry_transcripts.tsv   per-transcript read geometry
        read_geometry_summary.tsv       headline statistics, both ends
        read_geometry_tests.tsv         correlations and model coefficients
        read_geometry_strata.tsv        geometry binned by length and depth
        read_geometry_models.tsv        joint GLM of reach on length and depth
        read_geometry_shortfall.tsv     the shortfall prediction, tested per gene

Run:
  python3 read_geometry.py --bam DRS_our_ref_validation.bam \\
      --shortfall-genes Data/Shortfall/shortfall_genes.tsv --outdir Data/Shortfall
"""

import argparse
import os
import subprocess

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

AWK = r'''
/^@SQ/ {
    sn = ""; ln = 0
    for (i = 2; i <= NF; i++) {
        if (substr($i, 1, 3) == "SN:") sn = substr($i, 4)
        else if (substr($i, 1, 3) == "LN:") ln = substr($i, 4) + 0
    }
    L[sn] = ln
    next
}
/^@/ { next }
{
    cig = $6
    if (cig == "*") next
    num = 0; reflen = 0; clip5 = 0; clip3 = 0; firstop = 1; lastclip = 0
    n = length(cig)
    for (i = 1; i <= n; i++) {
        c = substr(cig, i, 1)
        if (c >= "0" && c <= "9") { num = num * 10 + (c - "0"); continue }
        if (c == "M" || c == "D" || c == "N" || c == "=" || c == "X") reflen += num
        if (c == "S" || c == "H") {
            if (firstop) clip5 = num
            lastclip = num
        } else {
            lastclip = 0
        }
        firstop = 0; num = 0
    }
    clip3 = lastclip
    ln = L[$3]
    if (ln == 0) next
    start_off = $4 - 1
    end_off = ln - ($4 + reflen - 1)
    if (end_off < 0) end_off = 0
    print $3 "\t" start_off "\t" end_off "\t" clip5 "\t" clip3 "\t" reflen "\t" ln
}
'''

NEAR = 10  # a read is "reaching" an end if it comes within this many nt


def q(a, p):
    return float(np.percentile(a, p)) if len(a) else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--shortfall-genes", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    OUT = args.outdir
    os.makedirs(OUT, exist_ok=True)

    print("PART 1. READ GEOMETRY, ALL PRIMARY SENSE ALIGNMENTS")
    print("samtools view -h -F 0x914: no secondary, supplementary, duplicate or reverse-strand alignment.")
    p1 = subprocess.Popen(["samtools", "view", "-h", "-F", "0x914", args.bam], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["awk", "-F", "\t", AWK], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    reads = pd.read_csv(p2.stdout, sep="\t", header=None,
                        names=["tx", "start_off", "end_off", "clip5", "clip3", "reflen", "tx_len"],
                        dtype={"tx": "category", "start_off": np.int32, "end_off": np.int32,
                               "clip5": np.int32, "clip3": np.int32, "reflen": np.int32, "tx_len": np.int32})
    p2.wait()
    if p2.returncode != 0:
        raise RuntimeError("samtools or awk failed")

    n_reads = len(reads)
    print("  alignments used %s   transcripts %s" % (f"{n_reads:,}", f"{reads['tx'].nunique():,}"))
    so, eo = reads["start_off"].to_numpy(), reads["end_off"].to_numpy()
    for name, a in (("5'", so), ("3'", eo)):
        print("  %s  median %.0f  <= %d nt: %.1f%%" % (name, q(a, 50), NEAR, 100 * float((a <= NEAR).mean())))

    print("\nPART 2. PER TRANSCRIPT GEOMETRY")
    g = reads.groupby("tx", observed=True)
    tx = pd.DataFrame({
        "n_reads": g.size(), "tx_len": g["tx_len"].first(), "med_reflen": g["reflen"].median(),
        "q05_start5": g["start_off"].quantile(0.05), "med_start5": g["start_off"].median(),
        "q05_end3": g["end_off"].quantile(0.05), "med_end3": g["end_off"].median(),
        "mean_clip5": g["clip5"].mean(), "mean_clip3": g["clip3"].mean(),
    })
    tx["frac_reach5"] = g["start_off"].apply(lambda s: float((s <= NEAR).mean()))
    tx["frac_reach3"] = g["end_off"].apply(lambda s: float((s <= NEAR).mean()))
    tx = tx.reset_index().rename(columns={"tx": "transcript"})
    tx["gene"] = tx["transcript"].astype(str).str.replace("_mRNA$", "", regex=True)
    tx["cov_frac"] = tx["med_reflen"] / tx["tx_len"]
    keep = tx[tx["n_reads"] >= 5].copy()
    print("  transcripts %s, with >=5 reads %s" % (f"{len(tx):,}", f"{len(keep):,}"))

    print("\nPART 3. DOES REACH DECAY WITH LENGTH, AND AT WHICH END")
    tests = []
    for lbl in ("frac_reach5", "frac_reach3"):
        r, p = stats.spearmanr(keep["tx_len"], keep[lbl])
        tests.append(("spearman", lbl, "tx_len", r, p, len(keep)))
        r2, p2v = stats.spearmanr(keep["n_reads"], keep[lbl])
        tests.append(("spearman", lbl, "n_reads", r2, p2v, len(keep)))
    r3, p3 = stats.spearmanr(keep["tx_len"], keep["med_start5"])
    tests.append(("spearman", "med_start5", "tx_len", r3, p3, len(keep)))
    r4, p4 = stats.spearmanr(keep["tx_len"], keep["med_end3"])
    tests.append(("spearman", "med_end3", "tx_len", r4, p4, len(keep)))
    for _, a, b, r, p, _n in tests:
        print("  %-14s vs %-10s rho=%.3f p=%.3g" % (a, b, r, p))

    strata = []
    keep["len_bin"] = pd.qcut(keep["tx_len"], 10, labels=False, duplicates="drop")
    for b, sub in keep.groupby("len_bin"):
        strata.append({"stratum": "tx_len", "bin": int(b) + 1, "lo": int(sub["tx_len"].min()),
                       "hi": int(sub["tx_len"].max()), "n": len(sub),
                       "mean_frac_reach5": round(float(sub["frac_reach5"].mean()), 4),
                       "mean_frac_reach3": round(float(sub["frac_reach3"].mean()), 4),
                       "median_med_start5": float(sub["med_start5"].median()),
                       "median_med_end3": float(sub["med_end3"].median())})
    keep["dep_bin"] = pd.qcut(keep["n_reads"], 10, labels=False, duplicates="drop")
    for b, sub in keep.groupby("dep_bin"):
        strata.append({"stratum": "n_reads", "bin": int(b) + 1, "lo": int(sub["n_reads"].min()),
                       "hi": int(sub["n_reads"].max()), "n": len(sub),
                       "mean_frac_reach5": round(float(sub["frac_reach5"].mean()), 4),
                       "mean_frac_reach3": round(float(sub["frac_reach3"].mean()), 4),
                       "median_med_start5": float(sub["med_start5"].median()),
                       "median_med_end3": float(sub["med_end3"].median())})

    print("\nPART 4. THE SHORTFALL PREDICTION, TESTED PER GENE")
    sg = pd.read_csv(args.shortfall_genes, sep="\t", dtype={"gene": str})
    m = keep.merge(sg, on="gene", how="inner")
    print("  genes joined to the shortfall table %s" % f"{len(m):,}")
    rows = []
    for end, flag, ntcol, obs in (("5'", "shortfall_5", "shortfall_nt_5", "q05_start5"),
                                  ("3'", "shortfall_3", "shortfall_nt_3", "q05_end3")):
        s = m[(m[flag] == 1) & m[ntcol].notna()].copy()
        if not len(s):
            continue
        r, p = stats.spearmanr(s[ntcol], s[obs])
        print("  %s end n=%d  predicted median=%.0f nt  observed median=%.0f nt  rho=%.3f p=%.3g"
              % (end, len(s), s[ntcol].median(), s[obs].median(), r, p))
        rows.append(("shortfall_edge", end, len(s), float(s[ntcol].median()), float(s[obs].median()), r, p))
        tests.append(("spearman", "q05_edge_%s" % end, ntcol, r, p, len(s)))
        ctl = m[(m[flag] == 0) & m[obs].notna()]
        u, pu = stats.mannwhitneyu(s[obs].dropna(), ctl[obs].dropna(), alternative="greater")
        n1, n2 = len(s[obs].dropna()), len(ctl[obs].dropna())
        rb = 2.0 * u / (n1 * n2) - 1.0
        tests.append(("mannwhitney", "q05_edge_%s" % end, "shortfall vs rest", rb, pu, n1 + n2))

    print("\nPART 5. THE CLASS shortfall.py CALLED UNAMBIGUOUS TRUNCATION")
    for end, flag, vcol, obs, reach, clip in (("5'", "shortfall_5", "verdict_5", "q05_start5", "frac_reach5", "mean_clip5"),
                                              ("3'", "shortfall_3", "verdict_3", "q05_end3", "frac_reach3", "mean_clip3")):
        s = m[(m[flag] == 1) & m[vcol].notna()].copy()
        if not len(s):
            continue
        for v, sub in s.groupby(vcol):
            print("  %s  %-16s n=%d reach=%.1f%% offset=%.0f clip=%.1f"
                  % (end, str(v), len(sub), 100 * sub[reach].mean(), sub[obs].median(), sub[clip].mean()))
            rows.append(("verdict_geometry", end, len(sub), str(v),
                        float(sub[reach].mean()), float(sub[obs].median()), float(sub[clip].mean())))
        vs = sorted(s[vcol].dropna().unique())
        if len(vs) >= 2:
            rr, p = stats.kruskal(*[s.loc[s[vcol] == v, reach].dropna() for v in vs])
            tests.append(("kruskal", "%s_reach_by_verdict" % end, "|".join(map(str, vs)), rr, p, len(s)))
        a = s.loc[s[vcol] == "both too short", clip].dropna()
        b = s.loc[s[vcol] != "both too short", clip].dropna()
        if len(a) and len(b):
            u, pu = stats.mannwhitneyu(a, b, alternative="greater")
            rb = 2.0 * u / (len(a) * len(b)) - 1.0
            tests.append(("mannwhitney", "%s_clip_both_too_short" % end, "rest of shortfall", rb, pu, len(a) + len(b)))

    print("\nPART 6. JOINT MODEL OF REACH")
    mods = []
    for end, col in (("5'", "frac_reach5"), ("3'", "frac_reach3")):
        f = keep.dropna(subset=[col, "tx_len", "n_reads"]).copy()
        X = sm.add_constant(pd.DataFrame({"log10_tx_len": np.log10(f["tx_len"].astype(float)),
                                          "log10_depth": np.log10(f["n_reads"].astype(float))}))
        res = sm.GLM(f[col].astype(float), X, family=sm.families.Binomial(),
                     freq_weights=f["n_reads"].astype(float)).fit()
        ci = res.conf_int()
        for term in ("log10_tx_len", "log10_depth"):
            b = res.params[term]
            lo, hi = ci.loc[term]
            print("  %s %-14s OR/decade=%.3f (%.3f-%.3f) p=%.3g" % (end, term, np.exp(b), np.exp(lo), np.exp(hi), res.pvalues[term]))
            mods.append({"model": "reach_%s" % end, "term": term, "n": len(f),
                        "or_per_decade": round(float(np.exp(b)), 4), "ci_lo": round(float(np.exp(lo)), 4),
                        "ci_hi": round(float(np.exp(hi)), 4), "p": float(res.pvalues[term])})

    cols = ["transcript", "gene", "n_reads", "tx_len", "med_reflen", "cov_frac",
            "q05_start5", "med_start5", "frac_reach5", "q05_end3", "med_end3", "frac_reach3",
            "mean_clip5", "mean_clip3"]
    tx[cols].round(4).to_csv(os.path.join(OUT, "read_geometry_transcripts.tsv"), sep="\t", index=False)
    summ = pd.DataFrame([
        {"end": "5'", "statistic": "reads", "value": n_reads},
        {"end": "5'", "statistic": "median_offset_nt", "value": q(so, 50)},
        {"end": "5'", "statistic": "pct_within_%dnt" % NEAR, "value": round(100.0 * float((so <= NEAR).mean()), 2)},
        {"end": "3'", "statistic": "reads", "value": n_reads},
        {"end": "3'", "statistic": "median_offset_nt", "value": q(eo, 50)},
        {"end": "3'", "statistic": "pct_within_%dnt" % NEAR, "value": round(100.0 * float((eo <= NEAR).mean()), 2)},
    ])
    summ.to_csv(os.path.join(OUT, "read_geometry_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(tests, columns=["test", "statistic", "against", "estimate", "p", "n"]).to_csv(
        os.path.join(OUT, "read_geometry_tests.tsv"), sep="\t", index=False)
    pd.DataFrame(rows, columns=["block", "end", "n", "group_or_predicted", "value_1", "value_2", "value_3"]).to_csv(
        os.path.join(OUT, "read_geometry_shortfall.tsv"), sep="\t", index=False)
    pd.DataFrame(strata).to_csv(os.path.join(OUT, "read_geometry_strata.tsv"), sep="\t", index=False)
    pd.DataFrame(mods).to_csv(os.path.join(OUT, "read_geometry_models.tsv"), sep="\t", index=False)
    print("\nwrote outputs to %s" % OUT)


if __name__ == "__main__":
    main()
