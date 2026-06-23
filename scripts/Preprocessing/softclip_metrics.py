#!/usr/bin/env python3
"""
softclip_metrics.py
-------------------
Computes per-transcript soft-clip statistics from a BAM file.
Designed for long-read (DRS / ONT) data aligned to a transcriptome reference.

For each transcript (reference sequence), reports:
  - n_reads          : number of primary alignments
  - mean_clip_total  : mean total soft-clip per read (5' + 3' combined)
  - mean_clip_5p     : mean soft-clip at the 5' end of the read
  - mean_clip_3p     : mean soft-clip at the 3' end of the read
  - frac_clipped     : fraction of reads with any soft-clip
  - mean_clip_frac   : mean (clip_bases / read_length) per read
  - asymmetry        : mean (clip_5p - clip_3p); positive = more 5' clipping
                       negative = more 3' clipping (UTR truncation signal)

Usage:
    python softclip_metrics.py --bam <file.bam> --out <output.tsv> [--min-reads 5]

Requirements:
    pysam (pip install pysam)
"""

import argparse
import sys
import pysam
import numpy as np
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True,
                        help="Input BAM file (must be sorted and indexed)")
    parser.add_argument("--out", required=True,
                        help="Output TSV file")
    parser.add_argument("--min-reads", type=int, default=5,
                        help="Minimum number of primary reads to report a transcript (default: 5)")
    parser.add_argument("--primary-only", action="store_true", default=True,
                        help="Use only primary alignments (default: True)")
    return parser.parse_args()


def get_softclip_lengths(read):
    """
    Parse CIGAR string and return (clip_5p, clip_3p) in bases.

    For a forward-strand read:
        - Leading S in CIGAR = 5' soft-clip
        - Trailing S in CIGAR = 3' soft-clip

    For a reverse-strand read the CIGAR is reported in read orientation,
    so leading S is still the read's 5' end, but maps to the transcript's
    3' end. We return clip relative to the READ orientation here;
    the caller handles strand flipping if needed.

    Returns (clip_5p, clip_3p) as integers.
    """
    if read.cigartuples is None:
        return 0, 0

    clip_5p = 0
    clip_3p = 0

    # pysam CIGAR op codes: 4 = soft clip (S)
    if read.cigartuples[0][0] == 4:
        clip_5p = read.cigartuples[0][1]
    if read.cigartuples[-1][0] == 4:
        clip_3p = read.cigartuples[-1][1]

    return clip_5p, clip_3p


def process_bam(bam_path, min_reads, primary_only):
    """
    Iterate over all reads in the BAM and accumulate per-transcript stats.
    Returns a dict: transcript_name -> dict of metric lists.
    """
    stats = defaultdict(lambda: {
        "clip_5p": [],
        "clip_3p": [],
        "clip_total": [],
        "clip_frac": [],
        "read_length": [],
    })

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            # Skip unmapped
            if read.is_unmapped:
                continue
            # Skip secondary and supplementary if primary_only
            if primary_only and (read.is_secondary or read.is_supplementary):
                continue

            transcript = read.reference_name
            clip_5p, clip_3p = get_softclip_lengths(read)

            # For reverse-strand reads, swap 5'/3' relative to transcript
            if read.is_reverse:
                clip_5p, clip_3p = clip_3p, clip_5p

            clip_total = clip_5p + clip_3p
            read_len = read.query_length if read.query_length else (read.infer_query_length() or 0)
            clip_frac = clip_total / read_len if read_len > 0 else 0.0

            stats[transcript]["clip_5p"].append(clip_5p)
            stats[transcript]["clip_3p"].append(clip_3p)
            stats[transcript]["clip_total"].append(clip_total)
            stats[transcript]["clip_frac"].append(clip_frac)
            stats[transcript]["read_length"].append(read_len)

    return stats


def summarise(stats, min_reads):
    """
    Convert per-read lists into per-transcript summary statistics.
    Returns list of dicts, one per transcript passing min_reads filter.
    """
    rows = []
    for transcript, d in stats.items():
        n = len(d["clip_total"])
        if n < min_reads:
            continue

        clip_total = np.array(d["clip_total"])
        clip_5p    = np.array(d["clip_5p"])
        clip_3p    = np.array(d["clip_3p"])
        clip_frac  = np.array(d["clip_frac"])

        rows.append({
            "transcript":       transcript,
            "n_reads":          n,
            "mean_read_length": round(np.mean(d["read_length"]), 1),
            "mean_clip_total":  round(np.mean(clip_total), 2),
            "median_clip_total":round(np.median(clip_total), 2),
            "mean_clip_5p":     round(np.mean(clip_5p), 2),
            "mean_clip_3p":     round(np.mean(clip_3p), 2),
            "frac_clipped":     round(np.mean(clip_total > 0), 4),
            "mean_clip_frac":   round(np.mean(clip_frac), 4),
            "asymmetry":        round(np.mean(clip_5p - clip_3p), 2),
        })

    # Sort by transcript name for reproducibility
    rows.sort(key=lambda r: r["transcript"])
    return rows


def write_output(rows, out_path):
    if not rows:
        print("WARNING: No transcripts passed the min-reads filter.", file=sys.stderr)
        return

    cols = list(rows[0].keys())
    with open(out_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    print(f"Written {len(rows)} transcripts to {out_path}", file=sys.stderr)


def main():
    args = parse_args()

    print(f"Processing: {args.bam}", file=sys.stderr)
    stats = process_bam(args.bam, args.min_reads, args.primary_only)

    print(f"Summarising {len(stats)} transcripts...", file=sys.stderr)
    rows = summarise(stats, args.min_reads)

    write_output(rows, args.out)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
