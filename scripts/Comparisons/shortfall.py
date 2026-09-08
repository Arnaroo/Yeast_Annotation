#!/usr/bin/env python3
"""
shortfall.py
==============
For genes where the DRS call is shorter than the prior reference, is the
transcript genuinely shorter, or is this DRS coverage bias at transcript
ends? TIF-seq (Pelechano et al., 2013) maps both transcript ends
directly, played no part in constructing the reference, and was never
used in the merge, so it can arbitrate. Underlies Supplementary Fig. S7
(plotted by plot_suppl_fig7.R, alongside read_geometry.py).

The mechanistic prediction is the second half of the answer. Direct RNA
sequencing is read 3' to 5'. If the 5' shortfall is a processivity
artefact it must grow with transcript length and shrink with read
depth, because both control the chance some read reaches the far end.
If the transcripts are genuinely shorter it should do neither. No such
mechanism exists at the 3' end, where sequencing starts.

COVARIATES
------------
Length is the spliced CDS length from the backbone GFF3, not the
transcript length -- regressing on transcript length would be circular,
since it contains the UTR call itself. Depth is the primary alignment
count against the ORF+1000nt reference from SRR32918129, an independent
DRS library (not one of the six construction libraries): a per-gene
coverage proxy for the protocol, not a count of the reads behind the
call.

Input:  Data/Annotation/{DRS_UTR_corrected.tsv, Nagalakshmi_UTR.csv,
        final_utr.tsv, Pelechano_UTR.csv, Pelechano_UTR.qc.csv,
        gene_models_backbone_with_UTR_slots.gff3}
        Data/Metrics/softclip/{DRS_orf_1000.softclip.tsv, DRS_our_ref.softclip.tsv}
Output: shortfall_summary.tsv       one row per end, the headline split
        shortfall_arbitration.tsv   TIF-seq verdict, three-way and two-way
        shortfall_genes.tsv         per-gene table, both ends, all covariates
        shortfall_models.tsv        logistic and rank model coefficients
        shortfall_strata.tsv        arbitration binned by length and by depth
        shortfall_merge_cost.tsv    what the merge rule costs against TIF-seq

Run:
  python3 shortfall.py --annotation-dir Data/Annotation \\
      --softclip-dir Data/Metrics/softclip --outdir Data/Shortfall
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

END_COL = {"5": "five", "3": "three"}
ENDS = ["5", "3"]


def log(msg=""):
    print(msg)


def fit_logit(frame, ycol, label, end):
    f = frame.dropna(subset=[ycol, "cds_len", "depth_orf1000"]).copy()
    f = f[f["depth_orf1000"] > 0]
    X = sm.add_constant(pd.DataFrame({
        "log10_cds_len": np.log10(f["cds_len"].astype(float)),
        "log10_depth": np.log10(f["depth_orf1000"].astype(float)),
    }))
    y = f[ycol].astype(float).to_numpy()
    res = sm.Logit(y, X).fit(disp=0)
    rows = []
    for term in ("log10_cds_len", "log10_depth"):
        b, se = res.params[term], res.bse[term]
        lo, hi = np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
        rows.append({
            "model": label, "end": end + "prime", "term": term, "n": len(f),
            "events": int(y.sum()), "odds_ratio_per_decade": round(float(np.exp(b)), 4),
            "ci_low": round(float(lo), 4), "ci_high": round(float(hi), 4),
            "p_value": float(res.pvalues[term]), "pseudo_r2": round(float(res.prsquared), 4),
        })
    log("  %s, %s' end   n=%d  events=%d (%.1f%%)  pseudoR2=%.4f"
        % (label, end, len(f), int(y.sum()), 100 * y.mean(), res.prsquared))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--softclip-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    ANN, SC, OUT = args.annotation_dir, args.softclip_dir, args.outdir
    os.makedirs(OUT, exist_ok=True)

    drs = pd.read_csv(os.path.join(ANN, "DRS_UTR_corrected.tsv"), sep="\t")
    drs.columns = ["gene", "five", "three"]
    nag = pd.read_csv(os.path.join(ANN, "Nagalakshmi_UTR.csv"))
    nag.columns = ["gene", "five", "three"]
    final = pd.read_csv(os.path.join(ANN, "final_utr.tsv"), sep="\t")
    final.columns = ["gene", "five", "three"]
    pel = pd.read_csv(os.path.join(ANN, "Pelechano_UTR.csv"))
    pel.columns = ["gene", "five", "three"]
    pel_qc = pd.read_csv(os.path.join(ANN, "Pelechano_UTR.qc.csv")).rename(columns={"Gene": "gene"})
    for df in (drs, nag, final, pel):
        for c in ("five", "three"):
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # -- CDS length, exogenous to every UTR quantity -------------------------
    id_re = re.compile(r"ID=([^;]+)")
    parent_re = re.compile(r"Parent=([^;]+)")
    gene_rows, cds_len = [], {}
    with open(os.path.join(ANN, "gene_models_backbone_with_UTR_slots.gff3")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[3].strip() == "NA" or f[4].strip() == "NA":
                continue
            start, end = int(f[3]), int(f[4])
            fid = id_re.search(f[8])
            if f[2] == "gene":
                if fid:
                    gene_rows.append((fid.group(1), f[0].strip(), f[6], start, end, end - start + 1))
            elif f[2] == "CDS" and fid and fid.group(1).endswith("_CDS"):
                m = parent_re.search(f[8])
                if m:
                    g = m.group(1).rsplit("_mRNA", 1)[0]
                    cds_len[g] = cds_len.get(g, 0) + (end - start + 1)
    genes = pd.DataFrame(gene_rows, columns=["gene", "chr", "strand", "gstart", "gend", "orf_span"])
    genes["cds_len"] = genes["gene"].map(cds_len)
    log("backbone genes: %d, with CDS length: %d" % (len(genes), int(genes["cds_len"].notna().sum())))

    sc = pd.read_csv(os.path.join(SC, "DRS_orf_1000.softclip.tsv"), sep="\t")
    sc = sc.rename(columns={"transcript": "gene", "n_reads": "depth_orf1000",
                            "mean_read_length": "mean_read_len"})[["gene", "depth_orf1000", "mean_read_len"]]
    scr = pd.read_csv(os.path.join(SC, "DRS_our_ref.softclip.tsv"), sep="\t")
    scr["gene"] = scr["transcript"].str.replace("_mRNA$", "", regex=True)
    scr = scr[["gene", "n_reads"]].rename(columns={"n_reads": "depth_ourref"})

    tab = genes[["gene", "chr", "strand", "orf_span", "cds_len"]].copy()
    for name, df in (("drs", drs), ("nag", nag), ("final", final), ("pel", pel)):
        tab = tab.merge(df.rename(columns={"five": name + "_5", "three": name + "_3"}), on="gene", how="left")
    tab = tab.merge(sc, on="gene", how="left").merge(scr, on="gene", how="left")
    tab = tab.merge(pel_qc[["gene", "mtif_support"]], on="gene", how="left")

    # ================= PART 1. The shortfall sets ==========================
    log("\nPART 1. WHERE DRS CALLS A SHORTER UTR THAN THE PRIOR REFERENCE")
    summary_rows = []
    for e in ENDS:
        d, n = tab["drs_" + e], tab["nag_" + e]
        both = d.notna() & n.notna()
        shorter, equal, longer = both & (d < n), both & (d == n), both & (d > n)
        mag = (n - d)[shorter]
        log("  %s' UTR: both=%d shorter=%d (%.1f%%) shortfall median=%.0f nt"
            % (e, both.sum(), shorter.sum(), 100 * shorter.sum() / both.sum(), mag.median()))
        summary_rows.append({
            "end": e + "prime", "n_both_called": int(both.sum()),
            "drs_shorter_n": int(shorter.sum()), "drs_shorter_pct": round(100 * shorter.sum() / both.sum(), 2),
            "identical_n": int(equal.sum()), "drs_longer_n": int(longer.sum()),
            "drs_longer_pct": round(100 * longer.sum() / both.sum(), 2),
            "shortfall_median_nt": float(mag.median()), "shortfall_q25_nt": float(mag.quantile(.25)),
            "shortfall_q75_nt": float(mag.quantile(.75)), "shortfall_max_nt": float(mag.max()),
            "drs_call_zero_n": int((d[shorter] == 0).sum()),
        })
        tab["shortfall_" + e] = np.where(both, shorter, np.nan)
        tab["shortfall_nt_" + e] = np.where(shorter, n - d, np.nan)
    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "shortfall_summary.tsv"), sep="\t", index=False)

    # ================= PART 2. TIF-seq as the arbiter ======================
    log("\nPART 2. WHICH CALL DOES THE INDEPENDENT BENCHMARK AGREE WITH")
    arb_rows = []
    for e in ENDS:
        d, n, p = tab["drs_" + e], tab["nag_" + e], tab["pel_" + e]
        sel = (tab["shortfall_" + e] == 1) & p.notna() & (p > 0)
        s = tab[sel]
        dd, nn, pp = s["drs_" + e], s["nag_" + e], s["pel_" + e]
        ed, en = (dd - pp).abs(), (nn - pp).abs()
        drs_closer, nag_closer, tied = int((ed < en).sum()), int((ed > en).sum()), int((ed == en).sum())
        below = int((pp <= dd).sum())
        between = int(((pp > dd) & (pp < nn)).sum())
        above = int((pp >= nn).sum())
        n_s = len(s)
        bt = stats.binomtest(drs_closer, drs_closer + nag_closer, 0.5)
        log("  %s' UTR n=%d  DRS closer %.1f%%  below %.1f%% between %.1f%% above %.1f%%"
            % (e, n_s, 100 * drs_closer / n_s, 100 * below / n_s, 100 * between / n_s, 100 * above / n_s))
        arb_rows.append({
            "end": e + "prime", "n": n_s, "subset": "drs_shorter",
            "median_drs_nt": float(dd.median()), "median_prior_nt": float(nn.median()),
            "median_tifseq_nt": float(pp.median()), "mae_drs_nt": float(ed.median()), "mae_prior_nt": float(en.median()),
            "drs_closer_n": drs_closer, "drs_closer_pct": round(100 * drs_closer / n_s, 2),
            "prior_closer_n": nag_closer, "prior_closer_pct": round(100 * nag_closer / n_s, 2),
            "tied_n": tied, "sign_test_p": float(bt.pvalue),
            "tifseq_below_both_n": below, "tifseq_below_both_pct": round(100 * below / n_s, 2),
            "tifseq_between_n": between, "tifseq_between_pct": round(100 * between / n_s, 2),
            "tifseq_above_both_n": above, "tifseq_above_both_pct": round(100 * above / n_s, 2),
        })
        tab.loc[sel, "verdict_" + e] = np.where(pp <= dd, "both too long",
                                                np.where(pp < nn, "DRS too short", "both too short"))

    log("\nCONTROL. THE SAME ARBITRATION WHERE DRS CALLS LONGER")
    for e in ENDS:
        d, n, p = tab["drs_" + e], tab["nag_" + e], tab["pel_" + e]
        sel = d.notna() & n.notna() & (d > n) & p.notna() & (p > 0)
        s = tab[sel]
        ed, en = (s["drs_" + e] - s["pel_" + e]).abs(), (s["nag_" + e] - s["pel_" + e]).abs()
        n_s = len(s)
        log("  %s' UTR (control) n=%d  DRS closer %.1f%%  prior closer %.1f%%"
            % (e, n_s, 100 * (ed < en).sum() / n_s, 100 * (ed > en).sum() / n_s))
        arb_rows.append({
            "end": e + "prime", "n": n_s, "subset": "drs_longer_control",
            "drs_closer_n": int((ed < en).sum()), "drs_closer_pct": round(100 * (ed < en).sum() / n_s, 2),
            "prior_closer_n": int((ed > en).sum()), "prior_closer_pct": round(100 * (ed > en).sum() / n_s, 2),
            "tied_n": int((ed == en).sum()), "mae_drs_nt": float(ed.median()), "mae_prior_nt": float(en.median()),
        })
    pd.DataFrame(arb_rows).to_csv(os.path.join(OUT, "shortfall_arbitration.tsv"), sep="\t", index=False)

    # ================= PART 3. The mechanistic test ========================
    log("\nPART 3. IS THE SHORTFALL A PROCESSIVITY ARTEFACT")
    model_rows = []
    for e in ENDS:
        model_rows += fit_logit(tab, "shortfall_" + e, "P(DRS shorter)", e)
    for e in ENDS:
        sub = tab[tab["verdict_" + e].notna()].copy()
        sub["truncating"] = (sub["verdict_" + e] == "both too short").astype(float)
        model_rows += fit_logit(sub, "truncating", "P(both too short | shortfall)", e)
    for e in ENDS:
        s = tab[tab["shortfall_" + e] == 1].dropna(subset=["shortfall_nt_" + e, "cds_len", "depth_orf1000"])
        s = s[s["depth_orf1000"] > 0]
        for cov in ("cds_len", "depth_orf1000"):
            r = stats.spearmanr(s["shortfall_nt_" + e], s[cov])
            model_rows.append({"model": "shortfall size", "end": e + "prime", "term": cov,
                               "n": len(s), "spearman_rho": round(float(r.statistic), 4), "p_value": float(r.pvalue)})
    pd.DataFrame(model_rows).to_csv(os.path.join(OUT, "shortfall_models.tsv"), sep="\t", index=False)

    # -- deciles, no model ----------------------------------------------------
    strat_rows = []
    for cov in ("cds_len", "depth_orf1000"):
        for e in ENDS:
            f = tab.dropna(subset=["shortfall_" + e, cov]).copy()
            f["bin"] = pd.qcut(f[cov], 10, labels=False, duplicates="drop")
            for b, g in f.groupby("bin"):
                v = g["verdict_" + e].dropna()
                trunc = (v == "both too short").mean() if len(v) else np.nan
                strat_rows.append({
                    "covariate": cov, "end": e + "prime", "decile": int(b) + 1,
                    "range_low": float(g[cov].min()), "range_high": float(g[cov].max()), "n": len(g),
                    "shortfall_pct": round(100 * float(g["shortfall_" + e].mean()), 2), "n_arbitrated": len(v),
                    "both_too_short_pct": round(100 * float(trunc), 2) if len(v) else np.nan,
                })
    pd.DataFrame(strat_rows).to_csv(os.path.join(OUT, "shortfall_strata.tsv"), sep="\t", index=False)

    # ================= PART 4. What the merge rule costs ===================
    log("\nPART 4. WHAT THE MERGE RULE DOES WITH THIS")
    cost_rows = []
    for e in ENDS:
        sel = tab["verdict_" + e].notna()
        s = tab[sel]
        kept_err = (s["nag_" + e] - s["pel_" + e]).abs()
        disc_err = (s["drs_" + e] - s["pel_" + e]).abs()
        delta = kept_err - disc_err
        w = stats.wilcoxon(kept_err, disc_err, alternative="two-sided", zero_method="wilcox")
        nz = delta[delta != 0]
        rb = float((np.sign(nz) * stats.rankdata(nz.abs())).sum() / (len(nz) * (len(nz) + 1) / 2))
        log("  %s' UTR n=%d  kept mae=%.0f discarded mae=%.0f penalty=%+.0f p=%.3g rb=%+.3f"
            % (e, len(s), kept_err.median(), disc_err.median(), delta.median(), w.pvalue, rb))
        cost_rows.append({
            "end": e + "prime", "n": len(s), "mae_kept_nt": float(kept_err.median()),
            "mae_discarded_nt": float(disc_err.median()), "median_penalty_nt": float(delta.median()),
            "wilcoxon_p": float(w.pvalue), "rank_biserial": round(rb, 4),
        })
    pd.DataFrame(cost_rows).to_csv(os.path.join(OUT, "shortfall_merge_cost.tsv"), sep="\t", index=False)

    keep = (["gene", "chr", "strand", "orf_span", "cds_len", "depth_orf1000", "depth_ourref",
             "mean_read_len", "mtif_support"]
            + [c + "_" + e for e in ENDS for c in ("drs", "nag", "final", "pel")]
            + ["shortfall_" + e for e in ENDS] + ["shortfall_nt_" + e for e in ENDS]
            + ["verdict_" + e for e in ENDS])
    tab[keep].to_csv(os.path.join(OUT, "shortfall_genes.tsv"), sep="\t", index=False)
    log("\nwrote shortfall_genes.tsv (%d genes)" % len(tab))


if __name__ == "__main__":
    main()
