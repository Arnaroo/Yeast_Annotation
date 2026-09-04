#!/usr/bin/env python3
"""
abundance_concordance.py
==========================
Tests whether changing the reference moves per-gene abundance estimates,
whether the genes that move are the ones whose extensions newly overlap
a neighbour, and whether that movement tracks the size of the boundary
change or looks like noise.

The estimator is a best-hit read count. The two references carry the
same transcripts under the same names and differ only in UTR length, so
the gene-level join is one-to-one and no read has to be apportioned
between competing isoforms of the same gene -- an EM-based
transcript-level quantifier is not required to ask this question.
Multi-mapping is reported rather than assumed away.

Counts are compared on two normalisations, which answer different
questions:
  CPM   library size only. "Does the gene capture more reads."
        The released reference is longer by construction, so part of
        any gain here is length.
  TPM   length normalised. "Does the estimate of abundance change once
        the extra length is accounted for." The stricter test, and the
        one a quantification pipeline would actually apply.
Both are reported; CPM alone would flatter the result.

Preparing the inputs
---------------------
For each of two technologies (direct RNA and short read) and two
references (prior, released), align the same reads with identical
settings and reduce to best-hit counts per transcript, e.g.

  minimap2 -ax map-ont  -k14 -N10 --secondary=no <ref.fa> <reads.fq.gz> \\
      | samtools view -b -F 0x904 - | samtools sort -o drs_new.bam
  samtools index drs_new.bam
  samtools idxstats drs_new.bam | awk '$3>0{print $1"\\t"$3}' > drs_new.counts.tsv

Repeat for drs_old, sr_new, sr_old against the prior reference. Also
needed: a two-column length table per reference (transcript, length in
nt), and Data/Overlap/neighbour_overlap_by_gene.tsv from
neighbour_overlap.py, for the "worst tier" (same-strand, into-ORF) test.

A length table is just each transcript's sequence length from the
corresponding FASTA already in Data/Annotation, e.g.

  seqkit fx2tab -nl Data/Annotation/Final_w.fa       > len_new.tsv
  seqkit fx2tab -nl Data/Annotation/Nagalakshmi_w.fa > len_old.tsv

Input:  {tech}_{new,old}.counts.tsv    two columns: transcript, count
        len_{new,old}.tsv              two columns: transcript, length
        Data/Annotation/final_utr.tsv, Nagalakshmi_UTR.csv
        Data/Overlap/neighbour_overlap_by_gene.tsv
Output: abundance_concordance.tsv   Spearman/Pearson, fold-change rates
        abundance_movers.tsv        worst-tier vs rest, Cliff's delta
        abundance_per_gene.tsv      full per-gene table, both technologies

Run:
  python3 abundance_concordance.py \\
      --counts-dir  work/counts \\
      --annotation-dir Data/Annotation \\
      --overlap-dir Data/Overlap \\
      --outdir      Data/Abundance
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--counts-dir", required=True,
                     help="directory holding {tech}_{new,old}.counts.tsv and len_{new,old}.tsv")
parser.add_argument("--annotation-dir", required=True)
parser.add_argument("--overlap-dir", required=True,
                     help="directory holding neighbour_overlap_by_gene.tsv (neighbour_overlap.py)")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

WORK, ANN, OVERLAP, OUT = args.counts_dir, args.annotation_dir, args.overlap_dir, args.outdir
os.makedirs(OUT, exist_ok=True)

FINAL_UTR = os.path.join(ANN, "final_utr.tsv")
NAGA_UTR = os.path.join(ANN, "Nagalakshmi_UTR.csv")
NEIGHBOUR = os.path.join(OVERLAP, "neighbour_overlap_by_gene.tsv")

# A 2-fold change on a log2 scale, fixed before looking at the data.
FC_THRESHOLD = 1.0
# Genes below this in either reference carry no usable estimate; the same
# floor used for the soft-clip analysis (five primary alignments).
MIN_COUNT = 5


def log(msg=""):
    print(msg)


def read_counts(label):
    path = os.path.join(WORK, f"{label}.counts.tsv")
    df = pd.read_csv(path, sep="\t", header=None, names=["transcript", "count"])
    df["gene"] = df["transcript"].str.replace(r"_mRNA$", "", regex=True)
    return df.set_index("gene")["count"]


def read_lengths(tag):
    path = os.path.join(WORK, f"len_{tag}.tsv")
    df = pd.read_csv(path, sep="\t", header=None, names=["transcript", "length"])
    df["gene"] = df["transcript"].str.replace(r"_mRNA$", "", regex=True)
    return df.set_index("gene")["length"]


def cliffs_delta(a, b):
    """Rank-based effect size, robust to heavy tails."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2.0 * u / (len(a) * len(b)) - 1.0


def main():
    log("abundance_concordance.py")
    log()

    lens = {"new": read_lengths("new"), "old": read_lengths("old")}
    universe = lens["new"].index
    log(f"transcripts in both references: {len(universe)}")

    tables = {}
    for tech in ("drs", "sr"):
        for tag in ("new", "old"):
            tables[f"{tech}_{tag}"] = read_counts(f"{tech}_{tag}").reindex(universe, fill_value=0)

    # ------------------------------------------- per-gene abundance concordance
    log("=" * 68)
    log("Per-gene abundance concordance")
    log("=" * 68)

    utr_final = pd.read_csv(FINAL_UTR, sep="\t").set_index("Gene")
    utr_naga = pd.read_csv(NAGA_UTR).set_index("Gene")
    ext = pd.DataFrame(index=universe)
    ext["len_new"] = lens["new"].reindex(universe)
    ext["len_old"] = lens["old"].reindex(universe)
    ext["length_gain_nt"] = ext["len_new"] - ext["len_old"]
    ext["extended"] = ext["length_gain_nt"] > 0

    neigh = pd.read_csv(NEIGHBOUR, sep="\t").set_index("Gene")

    per_tech = {}
    conc_rows = []
    for tech, techname in (("drs", "direct RNA"), ("sr", "short read")):
        c_new, c_old = tables[f"{tech}_new"], tables[f"{tech}_old"]
        df = pd.DataFrame({"count_new": c_new, "count_old": c_old})
        df["len_new"], df["len_old"] = ext["len_new"], ext["len_old"]
        df["length_gain_nt"], df["extended"] = ext["length_gain_nt"], ext["extended"]

        tot_new, tot_old = df["count_new"].sum(), df["count_old"].sum()
        df["cpm_new"] = 1e6 * df["count_new"] / tot_new
        df["cpm_old"] = 1e6 * df["count_old"] / tot_old
        rate_new = df["count_new"] / df["len_new"]
        rate_old = df["count_old"] / df["len_old"]
        df["tpm_new"] = 1e6 * rate_new / rate_new.sum()
        df["tpm_old"] = 1e6 * rate_old / rate_old.sum()

        keep = (df["count_new"] >= MIN_COUNT) & (df["count_old"] >= MIN_COUNT)
        df["quantifiable"] = keep
        sub = df[keep]

        for unit in ("cpm", "tpm"):
            lg = np.log2(sub[f"{unit}_new"]) - np.log2(sub[f"{unit}_old"])
            df.loc[keep, f"log2fc_{unit}"] = lg
            pear = stats.pearsonr(np.log2(sub[f"{unit}_new"]), np.log2(sub[f"{unit}_old"]))
            spear = stats.spearmanr(sub[f"{unit}_new"], sub[f"{unit}_old"])
            moved = lg.abs() > FC_THRESHOLD
            conc_rows.append({
                "technology": techname, "normalisation": unit.upper(),
                "n_quantifiable": int(keep.sum()), "pearson_log2": pear.statistic,
                "spearman": spear.statistic, "median_log2fc": lg.median(),
                "pct_within_1.2_fold": 100.0 * (lg.abs() < np.log2(1.2)).mean(),
                "n_moved_2fold": int(moved.sum()), "pct_moved_2fold": 100.0 * moved.mean(),
            })
        per_tech[tech] = df
        log(f"{techname}: {int(keep.sum()):,} of {len(df):,} genes reach "
            f"{MIN_COUNT} reads under both references")

    conc = pd.DataFrame(conc_rows)
    log()
    for _, r in conc.iterrows():
        log(f"  {r['technology']:11s} {r['normalisation']}  n={r['n_quantifiable']:>5,d}  "
            f"Spearman {r['spearman']:.4f}  within 1.2x {r['pct_within_1.2_fold']:5.1f}%  "
            f"moved 2x {r['n_moved_2fold']:>4,d} ({r['pct_moved_2fold']:.2f}%)")
    conc.to_csv(os.path.join(OUT, "abundance_concordance.tsv"), sep="\t", index=False)
    log()

    # -------------------------------------------- who moves, and is it the
    # -------------------------------------------- reviewer's newly-overlapping genes
    log("=" * 68)
    log("Newly-overlapping genes ('worst tier'), vs the rest")
    log("=" * 68)

    mover_rows = []
    for tech, techname in (("drs", "direct RNA"), ("sr", "short read")):
        df = per_tech[tech]
        sub = df[df["quantifiable"]].copy()
        sub["max_same_strand_orf_overlap_nt"] = (
            neigh["max_same_strand_orf_overlap_nt"].reindex(sub.index).fillna(0))
        sub["worst_tier"] = sub["max_same_strand_orf_overlap_nt"] > 0

        for unit in ("cpm", "tpm"):
            lg = sub[f"log2fc_{unit}"]
            moved = lg.abs() > FC_THRESHOLD
            a, b = lg[sub["worst_tier"]].abs(), lg[~sub["worst_tier"]].abs()
            delta_tier = cliffs_delta(a, b)
            mw = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else None
            mover_rows.append({
                "technology": techname, "normalisation": unit.upper(),
                "n_worst_tier": int(sub["worst_tier"].sum()),
                "median_abs_log2fc_worst_tier": a.median() if len(a) else np.nan,
                "median_abs_log2fc_rest": b.median() if len(b) else np.nan,
                "cliffs_delta_worst_tier": delta_tier,
                "mannwhitney_p_worst_tier": mw.pvalue if mw is not None else np.nan,
                "pct_moved_2fold_worst_tier": 100.0 * moved[sub["worst_tier"]].mean(),
                "pct_moved_2fold_rest": 100.0 * moved[~sub["worst_tier"]].mean(),
            })
        per_tech[tech] = sub

    movers = pd.DataFrame(mover_rows)
    for _, r in movers.iterrows():
        log(f"  {r['technology']:11s} {r['normalisation']}  n_worst_tier={r['n_worst_tier']:,}")
        log(f"    median |log2FC|  worst tier {r['median_abs_log2fc_worst_tier']:.4f}  "
            f"rest {r['median_abs_log2fc_rest']:.4f}  "
            f"Cliff delta {r['cliffs_delta_worst_tier']:+.3f}  "
            f"p {r['mannwhitney_p_worst_tier']:.3g}")
        log(f"    moved 2 fold      worst tier {r['pct_moved_2fold_worst_tier']:.2f}%  "
            f"rest {r['pct_moved_2fold_rest']:.2f}%")
    movers.to_csv(os.path.join(OUT, "abundance_movers.tsv"), sep="\t", index=False)
    log()

    keepcols = ["count_new", "count_old", "len_new", "len_old", "length_gain_nt",
                "extended", "quantifiable", "cpm_new", "cpm_old", "tpm_new", "tpm_old",
                "log2fc_cpm", "log2fc_tpm", "max_same_strand_orf_overlap_nt", "worst_tier"]
    out = []
    for tech, techname in (("drs", "direct RNA"), ("sr", "short read")):
        t = per_tech[tech].reindex(columns=keepcols).copy()
        t.insert(0, "technology", techname)
        out.append(t)
    pergene = pd.concat(out).reset_index().rename(columns={"index": "gene"})
    pergene.to_csv(os.path.join(OUT, "abundance_per_gene.tsv"), sep="\t", index=False)
    log(f"wrote abundance_per_gene.tsv, {len(pergene):,} rows")
    log("done")


if __name__ == "__main__":
    main()
