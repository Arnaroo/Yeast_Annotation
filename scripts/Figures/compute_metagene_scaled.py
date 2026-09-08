#!/usr/bin/env python3
"""
compute_metagene_scaled.py
=============================
Metagene coverage profiles for the our_ref BAM, split by whether the
final annotation extended that boundary by more than 20 nt relative to
Nagalakshmi, and scaled per gene so no small set of highly expressed
genes dominates the average. Underlies Supplementary Fig. S9 (plotted
by plot_suppl_fig9.py), a companion to Figure 4 restricted to a single
reference (our_ref) but adding the extension split and the scaled
statistic Figure 4 doesn't need.

Reviewer 2 asked for expression-normalised/per-gene-scaled profiles
alongside raw mean depth, and a version restricted to genes whose
boundaries were actually extended. Two defects had to be fixed first,
both also present in the raw Figure 4 profiles:

  A. Structural zeros. Padding the window with zero wherever the
     reference sequence runs out mixes "no signal" with "no sequence".
     Here a position is NaN where the sequence does not exist, the
     aggregation ignores NaN, and the number of contributing genes is
     written out at every position so a curve can be truncated where it
     stops being supported.
  B. Ratio instability. A mean of per-gene depth-ratios is dominated by
     the gene with the smallest denominator. Two guards: the scaled
     statistic only uses genes with an ORF-body depth of at least 5,
     and the median across genes (not the mean) is what gets plotted.

Anchors are the ORF start and ORF end, not the transcript start and
end -- the distinction Reviewer 2 comment 6 asked to be held throughout.

Input:  a BAM aligned to the our_ref transcriptome (reads filtered to
        primary, mapped alignments)
        Data/Annotation/final_utr.tsv, Nagalakshmi_UTR.csv
        a restrict-gene list, one systematic name per line (the shared
        gene set used throughout the validation, e.g. Data/Metrics/profiles/shared_genes.txt)
Output: metagene_scaled_profiles.tsv   one row per data_type x anchor x
                                       group(all/extended/unchanged) x
                                       statistic x position
        metagene_scaled_genes.tsv      per-gene body depth and extension flags

Run:
  python3 compute_metagene_scaled.py \\
      --bam DRS_our_ref.bam --data-type DRS \\
      --final-utr Data/Annotation/final_utr.tsv \\
      --nagalakshmi-utr Data/Annotation/Nagalakshmi_UTR.csv \\
      --restrict-genes Data/Metrics/profiles/shared_genes.txt \\
      --outdir Data/Metrics/profiles_scaled
  # repeat with --bam SR_merged_our_ref.bam --data-type SR --append
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import pysam

warnings.filterwarnings("ignore", r"All-NaN|Mean of empty|Degrees of freedom")

WINDOW = 1000
MIN_ORF = 200
THR = 20        # nt, same extension threshold as Fig. 2B
MIN_BODY = 5.0  # minimum ORF body depth to enter the scaled statistic
SUFFIX = "_mRNA"


def keep(r):
    """Read filter: primary, mapped alignments only."""
    return not r.is_secondary and not r.is_supplementary and not r.is_unmapped


def load_annotation(final_utr, nag_utr, restrict_path):
    fin = pd.read_csv(final_utr, sep="\t")
    nag = pd.read_csv(nag_utr)
    fin = fin.rename(columns={"Gene": "gene", "final_five_prime_utr": "fin5",
                              "final_three_prime_utr": "fin3"})
    nag = nag.rename(columns={"Gene": "gene", "five_prime_utr": "nag5",
                              "three_prime_utr": "nag3"})
    for c in ("fin5", "fin3"):
        fin[c] = pd.to_numeric(fin[c], errors="coerce")
    for c in ("nag5", "nag3"):
        nag[c] = pd.to_numeric(nag[c], errors="coerce")
    ann = fin.merge(nag, on="gene", how="left")
    with open(restrict_path) as fh:
        restrict = {l.strip() for l in fh if l.strip()}
    ann["in_restrict"] = ann["gene"].isin(restrict)
    ann["d5"] = ann["fin5"] - ann["nag5"]
    ann["d3"] = ann["fin3"] - ann["nag3"]
    ann["ext5"] = ann["d5"] > THR
    ann["ext3"] = ann["d3"] > THR
    return ann


def gene_profiles(bam_path, ann):
    lut = ann.set_index("gene")[["fin5", "fin3", "in_restrict"]].to_dict("index")
    genes, tss, tes, body, spans = [], [], [], [], []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        refs = list(zip(bam.references, bam.lengths))
        for ref, ref_len in refs:
            gene = ref[:-len(SUFFIX)] if ref.endswith(SUFFIX) else ref
            rec = lut.get(gene)
            if rec is None or not rec["in_restrict"]:
                continue
            u5, u3 = rec["fin5"], rec["fin3"]
            if pd.isna(u5) or pd.isna(u3):
                continue
            u5, u3 = int(u5), int(u3)
            orf_start, orf_end = u5, ref_len - u3
            orf_len = orf_end - orf_start
            if orf_len < MIN_ORF or orf_start < 0 or orf_end > ref_len:
                continue
            cov = np.array(bam.count_coverage(
                ref, quality_threshold=0, read_callback=keep)).sum(axis=0).astype(float)
            if len(cov) != ref_len:
                continue
            a = np.zeros(2 * WINDOW)
            up = min(WINDOW, orf_start)
            if up:
                a[WINDOW - up:WINDOW] = cov[orf_start - up:orf_start]
            dn = min(WINDOW, orf_len)
            if dn:
                a[WINDOW:WINDOW + dn] = cov[orf_start:orf_start + dn]
            b = np.zeros(2 * WINDOW)
            up = min(WINDOW, orf_len)
            if up:
                b[WINDOW - up:WINDOW] = cov[orf_end - up:orf_end]
            dn = min(WINDOW, ref_len - orf_end)
            if dn:
                b[WINDOW:WINDOW + dn] = cov[orf_end:orf_end + dn]
            genes.append(gene)
            tss.append(a)
            tes.append(b)
            body.append(cov[orf_start:orf_end].mean())
            spans.append((orf_start, orf_len, ref_len - orf_end))
    return (np.array(genes), np.vstack(tss), np.vstack(tes),
            np.array(body), np.array(spans))


def coverage_mask(spans, anchor):
    up = spans[:, 0] if anchor == "TSS" else spans[:, 1]
    dn = spans[:, 1] if anchor == "TSS" else spans[:, 2]
    pos = np.arange(-WINDOW, WINDOW)[None, :]
    return np.where(pos < 0, pos >= -up[:, None], pos < dn[:, None])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--data-type", required=True, choices=["DRS", "SR"])
    ap.add_argument("--final-utr", required=True)
    ap.add_argument("--nagalakshmi-utr", required=True)
    ap.add_argument("--restrict-genes", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--append", action="store_true",
                     help="append to existing outputs rather than overwrite "
                          "(run DRS first, then SR with --append)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ann = load_annotation(args.final_utr, args.nagalakshmi_utr, args.restrict_genes)
    dt = args.data_type
    print(f"{dt}  {args.bam}")
    g, tss, tes, body, spans = gene_profiles(args.bam, ann)
    print(f"  genes kept: {len(g)}")

    sub = ann.set_index("gene").loc[g]
    ext = {"TSS": sub["ext5"].to_numpy(), "TES": sub["ext3"].to_numpy()}
    ok = body >= MIN_BODY

    prof_rows = []
    pos = np.arange(-WINDOW, WINDOW)
    for anchor, mat in (("TSS", tss), ("TES", tes)):
        have = coverage_mask(spans, anchor)
        depth = np.where(have, mat, np.nan)
        ratio = np.where(have[ok], mat[ok] / body[ok][:, None], np.nan)

        for gname in ("all", "extended", "unchanged"):
            sel = (np.ones(len(g), bool) if gname == "all"
                   else ext[anchor] if gname == "extended" else ~ext[anchor])
            for stat in ("depth", "scaled_median"):
                if stat == "depth":
                    d = depth[sel]
                    cnt = np.isfinite(d).sum(axis=0)
                    val = np.nanmean(d, axis=0)
                    sem = np.nanstd(d, axis=0) / np.sqrt(np.maximum(cnt, 1))
                else:
                    d = ratio[sel[ok]]
                    cnt = np.isfinite(d).sum(axis=0)
                    val = np.nanmedian(d, axis=0)
                    lo, hi = np.nanpercentile(d, 25, axis=0), np.nanpercentile(d, 75, axis=0)
                    sem = (hi - lo) / 2.0
                prof_rows.append(pd.DataFrame({
                    "data_type": dt, "anchor": anchor, "group": gname,
                    "statistic": stat, "position": pos, "value": val, "spread": sem,
                    "n_genes": int(sel.sum() if stat == "depth" else sel[ok].sum()),
                    "n_at_position": cnt,
                }))
        print(f"  {anchor}  extended {int(ext[anchor].sum())}  "
              f"unchanged {int((~ext[anchor]).sum())}  scaled set {int(ok.sum())}")

    prof = pd.concat(prof_rows, ignore_index=True)
    gt = pd.DataFrame({"gene": g, "body_mean_depth": body,
                       "d5": sub["d5"].to_numpy(), "d3": sub["d3"].to_numpy(),
                       "ext5": ext["TSS"], "ext3": ext["TES"]}).assign(data_type=dt)

    prof_path = os.path.join(args.outdir, "metagene_scaled_profiles.tsv")
    gene_path = os.path.join(args.outdir, "metagene_scaled_genes.tsv")
    if args.append and os.path.exists(prof_path):
        prof = pd.concat([pd.read_csv(prof_path, sep="\t"), prof], ignore_index=True)
        gt = pd.concat([pd.read_csv(gene_path, sep="\t"), gt], ignore_index=True)
    prof.to_csv(prof_path, sep="\t", index=False)
    gt.to_csv(gene_path, sep="\t", index=False)
    print(f"wrote {prof_path} ({len(prof)} rows)")
    print(f"wrote {gene_path} ({len(gt)} rows)")


if __name__ == "__main__":
    main()
