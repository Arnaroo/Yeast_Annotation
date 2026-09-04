#!/usr/bin/env python3
"""
read_migration_analyse.py
===========================
Turns the per-technology transition tables from read_migration.sh into
the read-level answer to whether an extended UTR reaching into a
neighbour takes reads that belong to the neighbour.

abundance_concordance.py shows abundance is nearly unchanged in
aggregate, which is reassuring but not the same claim: two genes can
swap a thousand reads and still sit on the diagonal if the swap is
roughly symmetric. Following individual reads is the only way to see a
directional theft directly.

Categories
----------
  both_unmapped   unmapped under both references
  unchanged       same gene under both
  gained          unmapped under the prior reference, assigned under
                  the released one -- the intended effect of a longer
                  UTR, not a cost
  lost            assigned under the prior reference, unmapped under
                  the released one
  migrated        assigned to two different genes

Migration is then split three ways, which is where the argument is won
or lost:
  attributable        the gene pair is one of the newly overlapping
                       pairs from neighbour_overlap.py -- the only
                       migrations this work could have caused
  attributable_into_orf  same, restricted to pairs where the overlap
                       reaches into the partner's ORF body on the same
                       strand -- where a read is genuinely ambiguous
  other                everything else: multi-mapping jitter between
                       similar sequences, which would occur between any
                       two references and is not a property of this one

The attributable count is an upper bound on harm, not an estimate: a
read landing on a newly overlapping pair is not proof the new
assignment is wrong.

The RDN37 rDNA repeat (YLR154W-A..H / YLR154C-G/H) is a known
multi-mapping sink for every yeast reference, present in the assembly
at reduced copy number, and is reported separately from the rest of the
attributable migration for that reason.

Input:  <transitions-dir>/{drs,sr}.transitions.tsv  (read_migration.sh)
        <overlap-dir>/neighbour_overlap_pairs.tsv    (neighbour_overlap.py)
Output: read_migration_summary.tsv   one row per technology x measure
        read_migration_pairs.tsv     top attributable gene-pair transitions

Run:
  python3 read_migration_analyse.py \
      --transitions-dir work/migration --overlap-dir Data/Overlap --outdir Data/Abundance
"""

import argparse
import os

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--transitions-dir", required=True)
parser.add_argument("--overlap-dir", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

TRANS, OVERLAP, OUT = args.transitions_dir, args.overlap_dir, args.outdir
os.makedirs(OUT, exist_ok=True)

UNMAPPED = "UNMAPPED"

# YLR154W-A to -H sit inside the RDN37 rDNA repeat on chromosome XII,
# present in the assembly at reduced copy number and a known
# multi-mapping sink for every yeast reference, ours and the prior one
# alike. Separated out rather than removed: the total is still reported,
# but a total dominated by one repeat locus means something different
# from the same total spread evenly.
RDNA_LOCUS = {f"YLR154W-{s}" for s in "ABCDEFGH"} | {"YLR154C-G", "YLR154C-H"}


def log(msg=""):
    print(msg)


def gene_of(name):
    if name == UNMAPPED:
        return UNMAPPED
    return name[:-5] if name.endswith("_mRNA") else name


def main():
    pairs = pd.read_csv(os.path.join(OVERLAP, "neighbour_overlap_pairs.tsv"), sep="\t")
    new_overlap = {frozenset((a, b)) for a, b in zip(pairs["gene_a"], pairs["gene_b"])}
    log(f"newly overlapping pairs: {len(new_overlap)}")

    worst = pairs[(pairs["same_strand"].astype(str) == "True") & (pairs["overlap_into_orf_nt"] > 0)]
    worst_pairs = {frozenset((a, b)) for a, b in zip(worst["gene_a"], worst["gene_b"])}
    log(f"  of which same strand and reaching an ORF body: {len(worst_pairs)}")
    log()

    summary_rows, pair_rows = [], []

    for tech, label in (("drs", "direct RNA"), ("sr", "short read")):
        path = os.path.join(TRANS, f"{tech}.transitions.tsv")
        if not os.path.exists(path):
            log(f"SKIP {tech}, {path} absent")
            continue

        t = pd.read_csv(path, sep="\t")
        total = int(t["n_reads"].sum())
        log(f"{label}  {len(t)} distinct transitions, {total:,} reads")

        t["gene_old"] = t["old"].map(gene_of)
        t["gene_new"] = t["new"].map(gene_of)
        old_un, new_un = t["gene_old"].eq(UNMAPPED), t["gene_new"].eq(UNMAPPED)
        same = t["gene_old"].eq(t["gene_new"])

        cat = pd.Series("migrated", index=t.index)
        cat[old_un & new_un] = "both_unmapped"
        cat[~old_un & ~new_un & same] = "unchanged"
        cat[old_un & ~new_un] = "gained"
        cat[~old_un & new_un] = "lost"
        t["category"] = cat

        counts = t.groupby("category")["n_reads"].sum()
        for c in ("both_unmapped", "unchanged", "gained", "lost", "migrated"):
            n = int(counts.get(c, 0))
            log(f"  {c:<14} {n:>12,}  {100.0 * n / total:6.3f} %")
            summary_rows.append(dict(technology=label, measure=c, n_reads=n,
                                      pct_of_reads=round(100.0 * n / total, 4)))

        mig = t[t["category"] == "migrated"].copy()
        n_mig = int(mig["n_reads"].sum())
        if n_mig == 0:
            log("  no migration, nothing to split")
            log()
            continue

        key = [frozenset((a, b)) for a, b in zip(mig["gene_old"], mig["gene_new"])]
        mig["attributable"] = [k in new_overlap for k in key]
        mig["worst_tier"] = [k in worst_pairs for k in key]

        n_attr = int(mig.loc[mig["attributable"], "n_reads"].sum())
        n_worst = int(mig.loc[mig["worst_tier"], "n_reads"].sum())
        n_other = n_mig - n_attr

        log(f"  migration split, denominator is all {total:,} reads")
        for name, n in (("migrated_attributable", n_attr),
                        ("migrated_attributable_into_orf", n_worst),
                        ("migrated_other", n_other)):
            log(f"    {name:<32} {n:>10,}  {100.0 * n / total:6.4f} %")
            summary_rows.append(dict(technology=label, measure=name, n_reads=n,
                                      pct_of_reads=round(100.0 * n / total, 4)))

        att_all = mig[mig["attributable"]]
        at_rdna = att_all["gene_old"].isin(RDNA_LOCUS) | att_all["gene_new"].isin(RDNA_LOCUS)
        n_rdna = int(att_all.loc[at_rdna, "n_reads"].sum())
        n_else = int(att_all.loc[~at_rdna, "n_reads"].sum())
        elsewhere = att_all[~at_rdna]
        n_genes_else = len(set(elsewhere["gene_old"]) | set(elsewhere["gene_new"]))
        for name, n in (("migrated_attributable_rdna_repeat", n_rdna),
                        ("migrated_attributable_elsewhere", n_else)):
            log(f"    {name:<32} {n:>10,}  {100.0 * n / total:6.4f} %")
            summary_rows.append(dict(technology=label, measure=name, n_reads=n,
                                      pct_of_reads=round(100.0 * n / total, 4)))
        log(f"    {'genes involved elsewhere':<32} {n_genes_else:>10,}")

        att = mig[mig["attributable"]].sort_values("n_reads", ascending=False)
        for _, r in att.head(200).iterrows():
            pair_rows.append(dict(technology=label, gene_old=r["gene_old"], gene_new=r["gene_new"],
                                   n_reads=int(r["n_reads"]),
                                   pct_of_reads=round(100.0 * r["n_reads"] / total, 6),
                                   into_orf_body=bool(r["worst_tier"])))
        if len(att):
            top = att.iloc[0]
            log(f"    worst single pair: {top['gene_old']} to {top['gene_new']}, "
                f"{int(top['n_reads']):,} reads")
        log()

    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "read_migration_summary.tsv"),
                                       sep="\t", index=False)
    pd.DataFrame(pair_rows).to_csv(os.path.join(OUT, "read_migration_pairs.tsv"),
                                    sep="\t", index=False)
    log("wrote read_migration_summary.tsv")
    log("wrote read_migration_pairs.tsv")


if __name__ == "__main__":
    main()
