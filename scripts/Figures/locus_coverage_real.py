#!/usr/bin/env python3
"""
locus_coverage_real.py
=========================
Real per-position coverage for the nine Figure 1B / Supplementary
Fig. S3 windows defined by locus_windows.py, summed over the four
fraction-level construction conditions (NS_PS, NS_RDT, S10_PS, S10_RDT
-- not the NS/S10 pooled tables, which are themselves those four summed
and would double count).

ANTISENSE TRACK
------------------
Each condition table is single-stranded by construction (one gene's
own window, aligned against a sense-oriented reference for that gene
alone), so a gene's own table carries no antisense signal. What a
locus plot shows as antisense is a NEIGHBOURING gene's own sense
signal, falling inside this gene's window because the two ORFs are
close: for every opposite-strand gene whose own ORF+-1000nt window
overlaps this panel (read from locus_features.tsv), that gene's own
real coverage is pulled, mapped into genomic coordinates by its own
span and strand, and summed into the antisense track.

WINDOW-POSITION COLUMN, AND WHY IT MATTERS FOR MINUS-STRAND GENES
---------------------------------------------------------------------
Each condition table indexes position 1..n in the gene's own 5'->3'
transcript direction (matching the segmentation pipeline exactly), not
ascending genomic order. This script keeps a `window_position` column
(1..n in the gene's own transcript direction) alongside the genomic
`position` column, and locus_segmentation.py must segment on
window_position, not on genomic order -- for a minus-strand gene those
run in opposite directions, and segmenting on genomic order silently
swaps the five/three sides.

Input:  locus_windows.tsv, locus_features.tsv (locus_windows.py)
        {condition}_combined.tsv for NS_PS, NS_RDT, S10_PS, S10_RDT
        (ReferenceConstruction/01_aggregate_data.R), columns gene,
        position, transcript_np
Output: locus_coverage.tsv    per position, per strand depth, plus
                              window_position
        locus_windows.tsv     rewritten in place with a coverage-source column

Run:
  python3 locus_coverage_real.py --indir Data/LocusTracks \\
      --combined-dir Data/CombinedProfiles --outdir Data/LocusTracks
"""

import argparse
import os
import subprocess

import numpy as np
import pandas as pd

CONDITIONS = ["NS_PS", "NS_RDT", "S10_PS", "S10_RDT"]
PAD = 1000


def gene_window_table(gene, combined_dir):
    """Sum transcript_np for `gene` across the four condition tables."""
    total = None
    for cond in CONDITIONS:
        path = os.path.join(combined_dir, f"{cond}_combined.tsv")
        # grep avoids loading a many-million-row table just to keep one gene's
        # rows; anchored at line start with a real tab so a gene name that is
        # a suffix of another cannot match the wrong rows.
        with open(path) as fh:
            out = subprocess.run(["grep", "-P", f"^{gene}\t"], stdin=fh, capture_output=True, text=True)
        if not out.stdout:
            continue
        rows = [line.split("\t") for line in out.stdout.splitlines()]
        pos = np.array([int(r[1]) for r in rows], dtype=np.int64)
        val = np.array([int(r[2]) for r in rows], dtype=np.int64)
        order = np.argsort(pos)
        pos, val = pos[order], val[order]
        if total is None:
            total = np.zeros(pos.max(), dtype=np.int64)
        if pos.max() > len(total):
            total = np.pad(total, (0, pos.max() - len(total)))
        total[pos - 1] += val
    return total


def to_genomic(depth_by_pos, g_lo, g_hi, strand):
    n = len(depth_by_pos)
    out = {}
    for p in range(1, n + 1):
        g = (g_lo - PAD + (p - 1)) if strand == "+" else (g_hi + PAD - (p - 1))
        out[g] = out.get(g, 0) + int(depth_by_pos[p - 1])
    return out


def to_genomic_with_window_pos(depth_by_pos, g_lo, g_hi, strand):
    n = len(depth_by_pos)
    depth, wpos = {}, {}
    for p in range(1, n + 1):
        g = (g_lo - PAD + (p - 1)) if strand == "+" else (g_hi + PAD - (p - 1))
        depth[g] = depth.get(g, 0) + int(depth_by_pos[p - 1])
        wpos[g] = p
    return depth, wpos


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True, help="directory holding locus_windows.tsv, locus_features.tsv")
    ap.add_argument("--combined-dir", required=True,
                     help="directory holding {condition}_combined.tsv for the four conditions")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    win = pd.read_csv(os.path.join(args.indir, "locus_windows.tsv"), sep="\t")
    feat = pd.read_csv(os.path.join(args.indir, "locus_features.tsv"), sep="\t")
    gene_cache = {}

    def cached_window(gene):
        if gene not in gene_cache:
            gene_cache[gene] = gene_window_table(gene, args.combined_dir)
        return gene_cache[gene]

    cov_rows = []
    for _, w in win.iterrows():
        panel, gene, chrom, strand = w["panel"], w["gene"], w["chrom"], w["strand"]
        lo, hi = int(w["window_start"]), int(w["window_end"])
        g_lo, g_hi = int(w["orf_start"]), int(w["orf_end"])

        sense_win = cached_window(gene)
        if sense_win is None:
            print(f"WARNING {panel} {gene}: no real coverage found, sense track is all zero")
            sense_g, wpos_g = {}, {}
        else:
            sense_g, wpos_g = to_genomic_with_window_pos(sense_win, g_lo, g_hi, strand)

        neigh = feat[(feat.panel == panel) & (feat.feature_class == "gene") &
                     (feat.source == "final") & (feat.strand != strand)]
        anti_g, neigh_used, neigh_missing = {}, [], []
        for _, nrow in neigh.iterrows():
            ng, ns, ng_lo, ng_hi = nrow["feature_gene"], nrow["strand"], int(nrow["start"]), int(nrow["end"])
            nwin = cached_window(ng)
            if nwin is None:
                neigh_missing.append(ng)
                continue
            neigh_used.append(ng)
            for g, v in to_genomic(nwin, ng_lo, ng_hi, ns).items():
                if lo <= g <= hi:
                    anti_g[g] = anti_g.get(g, 0) + v

        pos = np.arange(lo, hi + 1, dtype=np.int64)
        sense = np.array([sense_g.get(int(p), 0) for p in pos], dtype=np.int64)
        anti = np.array([anti_g.get(int(p), 0) for p in pos], dtype=np.int64)
        wpos = np.array([wpos_g.get(int(p), -1) for p in pos], dtype=np.int64)

        cov_rows.append(pd.DataFrame({"panel": panel, "gene": gene, "chrom": chrom, "position": pos,
                                      "window_position": wpos, "depth_sense": sense, "depth_antisense": anti}))
        print(f"{panel:6s} {gene:12s} {chrom}:{lo}-{hi} ({strand})  max sense {sense.max():6d}  "
              f"max antisense {anti.max():6d}  antisense from {neigh_used or 'none'}"
              + (f"  MISSING coverage for {neigh_missing}" if neigh_missing else ""))

    pd.concat(cov_rows, ignore_index=True).to_csv(os.path.join(args.outdir, "locus_coverage.tsv"), sep="\t", index=False)
    win2 = win.copy()
    win2["coverage_source"] = f"real ({'+'.join(CONDITIONS)}, summed per-condition combined tables)"
    win2.to_csv(os.path.join(args.outdir, "locus_windows.tsv"), sep="\t", index=False)
    print("\nwrote locus_coverage.tsv")
    print("wrote locus_windows.tsv (coverage-source column added)")


if __name__ == "__main__":
    main()
