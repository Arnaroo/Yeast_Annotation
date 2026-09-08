#!/usr/bin/env python3
"""
locus_windows.py
===================
Defines the genomic windows and annotation tracks for Figure 1B and
Supplementary Fig. S3: each window is the gene's ORF span padded by
1000 nt each side (the same window 02_transcripts_segmentation.R hands
to the segmenter), plus every CDS/UTR/intron feature from both the
final and Nagalakshmi GFF3s that overlaps it, so neighbouring genes
appear on the CDS track.

Coverage itself is NOT computed here -- see locus_coverage_real.py,
which reuses this script's window/feature tables and reads the real
per-condition coverage.

Input:  Data/Annotation/Final_UTRs.gff3, Nagalakshmi_UTRs.gff3
Output: locus_features.tsv   CDS, UTR and intron spans, both sources
        locus_windows.tsv    one row per panel, window definition (no
                              coverage-dependent columns yet)

Run:
  python3 locus_windows.py --annotation-dir Data/Annotation --outdir Data/LocusTracks
"""

import argparse
import os
from collections import defaultdict

import pandas as pd

PAD = 1000  # 02_transcripts_segmentation.R's borne_inf/borne_sup

# Figure 1 panel B locus, then the eight Supplementary Fig. S3 loci in
# the order the caption lists them (chosen to be difficult, not
# representative: a close neighbour, or an intron-containing gene).
PANELS = [
    ("Fig1B", "YBL091C"),
    ("S1A", "YFR008W"), ("S1B", "YOR110W"), ("S1C", "YBR084W"), ("S1D", "YCL045C"),
    ("S1E", "YCR015C"), ("S1F", "YDL167C"), ("S1G", "YDR376W"), ("S1H", "YGL103W"),
]


def parse_gff(path):
    """
    Return {gene: [(chrom, feature_class, start, end, strand), ...]}.

    The GFF3 files write UTRs and introns with a CDS feature type and
    encode the real class in the ID attribute, so class is taken from
    the ID. Genes with no UTR call at an end carry the literal string
    NA in the coordinate columns rather than being omitted; those
    records are skipped and the count reported.
    """
    feats = defaultdict(list)
    n_na = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = p[:9]
            if start.strip() == "NA" or end.strip() == "NA":
                n_na += 1
                continue
            start, end = int(start.strip()), int(end.strip())
            fid = ""
            for kv in attrs.split(";"):
                if kv.startswith("ID="):
                    fid = kv[3:]
            if fid.endswith("_five_prime_UTR"):
                gene, cls = fid[:-15], "five_prime_UTR"
            elif fid.endswith("_three_prime_UTR"):
                gene, cls = fid[:-16], "three_prime_UTR"
            elif fid.endswith("_CDS"):
                gene, cls = fid[:-4], "CDS"
            elif fid.endswith("_intron"):
                gene, cls = fid[:-7], "intron"
            elif ftype == "gene":
                gene, cls = fid, "gene"
            else:
                continue
            feats[gene].append((chrom, cls, start, end, strand))
    print(f"  {os.path.basename(path)}: skipped {n_na} records with NA coordinates, {len(feats)} genes")
    return feats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    final = parse_gff(os.path.join(args.annotation_dir, "Final_UTRs.gff3"))
    naga = parse_gff(os.path.join(args.annotation_dir, "Nagalakshmi_UTRs.gff3"))

    feat_rows, win_rows = [], []
    for panel, gene in PANELS:
        orf = [r for r in final[gene] if r[1] == "gene"]
        if not orf:
            raise SystemExit(f"{gene}: no gene record in Final_UTRs.gff3")
        chrom, _, g_lo, g_hi, strand = orf[0]
        lo, hi = max(1, g_lo - PAD), g_hi + PAD  # clamped properly once chrom lengths are known

        n_feat = 0
        for src, table in (("final", final), ("nagalakshmi", naga)):
            for other, rl in table.items():
                for c, cls, s, e, st in rl:
                    if c != chrom or e < lo or s > hi:
                        continue
                    feat_rows.append({"panel": panel, "window_gene": gene, "source": src,
                                      "feature_gene": other, "feature_class": cls, "chrom": c,
                                      "start": s, "end": e, "strand": st, "is_window_gene": int(other == gene)})
                    n_feat += 1

        win_rows.append({"panel": panel, "gene": gene, "chrom": chrom, "strand": strand,
                         "orf_start": g_lo, "orf_end": g_hi, "window_start": lo, "window_end": hi,
                         "window_len": hi - lo + 1, "n_features_drawn": n_feat})
        print(f"{panel:6s} {gene:12s} {chrom}:{lo}-{hi} ({strand})  features {n_feat}")

    pd.DataFrame(feat_rows).to_csv(os.path.join(args.outdir, "locus_features.tsv"), sep="\t", index=False)
    pd.DataFrame(win_rows).to_csv(os.path.join(args.outdir, "locus_windows.tsv"), sep="\t", index=False)
    print("\nwrote locus_features.tsv")
    print("wrote locus_windows.tsv")


if __name__ == "__main__":
    main()
