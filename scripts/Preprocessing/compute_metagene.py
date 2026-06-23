#!/usr/bin/env python3
"""
compute_metagene.py  — v4
--------------------------
Computes metagene coverage profiles anchored at the true ORF start and
ORF end, from transcriptome BAM files. Coverage is shown everywhere the
reference sequence exists — i.e. ORF-only has zero coverage outside the
ORF (no sequence there), ORF±1000 has real coverage in its flanks, and
our_ref / Nagalakshmi show fading UTR coverage that goes to zero once
the window exceeds the annotated UTR length.

GENE RESTRICTION (--restrict_to):
  Pass a file of gene names (one per line, no suffix) to restrict ALL
  four references to the exact same transcript set. This is required
  for a fair comparison: if Nagalakshmi only has UTR data for a subset
  of genes, comparing its metagene against an unrestricted orf_only/
  orf_1000/our_ref metagene compares different gene universes and the
  result is not interpretable.

USAGE — orf_only (no UTRs, contig names have NO suffix):
    python compute_metagene.py \\
        --bam DRS_orf_only.bam --outdir profiles/ \\
        --restrict_to shared_genes.txt

USAGE — orf_1000 (contig names have NO suffix):
    python compute_metagene.py \\
        --bam DRS_orf_1000.bam --fixed_utr5 1000 --fixed_utr3 1000 \\
        --outdir profiles/ --restrict_to shared_genes.txt

USAGE — our_ref (contig names end in _mRNA):
    python compute_metagene.py \\
        --bam DRS_our_ref.bam \\
        --utr_file final_utr.tsv \\
        --utr5_col final_five_prime_utr --utr3_col final_three_prime_utr \\
        --name_suffix _mRNA \\
        --outdir profiles/ --restrict_to shared_genes.txt

USAGE — nagalakshmi (contig names end in _mRNA):
    python compute_metagene.py \\
        --bam DRS_nagalakshmi.bam \\
        --utr_file Nagalakshmi_UTR.csv \\
        --utr5_col five_prime_utr --utr3_col three_prime_utr \\
        --name_suffix _mRNA \\
        --outdir profiles/ --restrict_to shared_genes.txt

Requirements: pysam, numpy, pandas
"""

import argparse
import sys
import numpy as np
import pandas as pd
import pysam
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam",          required=True)
    p.add_argument("--outdir",       required=True)
    p.add_argument("--window",       type=int, default=1000)
    p.add_argument("--min_orf",      type=int, default=200)
    p.add_argument("--utr_file",     default=None,
                   help="TSV/CSV with per-gene UTR lengths (our_ref/nagalakshmi)")
    p.add_argument("--utr5_col",     default=None)
    p.add_argument("--utr3_col",     default=None)
    p.add_argument("--gene_col",     default=None,
                   help="Gene name column (default: first column)")
    p.add_argument("--fixed_utr5",   type=int, default=0,
                   help="Fixed 5'UTR length, used when --utr_file is not given")
    p.add_argument("--fixed_utr3",   type=int, default=0,
                   help="Fixed 3'UTR length, used when --utr_file is not given")
    p.add_argument("--name_suffix",  default="",
                   help="Suffix on BAM contig names to strip before gene "
                        "lookup, e.g. '_mRNA'. Leave empty if contigs are "
                        "already bare gene names (orf_only, orf_1000).")
    p.add_argument("--restrict_to",  default=None,
                   help="File of gene names (one per line, bare, no suffix) "
                        "to restrict ALL references to the same gene set.")
    return p.parse_args()


def strip_suffix(name, suffix):
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    return name


def load_utr_table(path, utr5_col, utr3_col, gene_col=None):
    """Return dict: gene_name -> (utr5, utr3) or None if NA/missing."""
    sep = "," if path.endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep)
    if gene_col is None:
        gene_col = df.columns[0]
    out = {}
    for _, row in df.iterrows():
        gene = str(row[gene_col])
        try:
            u5 = int(float(row[utr5_col]))
            u3 = int(float(row[utr3_col]))
        except (ValueError, TypeError):
            out[gene] = None
            continue
        out[gene] = (u5, u3)
    return out


def load_restrict_set(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def compute_profiles(bam_path, utr_lookup, fixed_utr5, fixed_utr3,
                     name_suffix, window, min_orf, restrict_genes):
    """
    Iterate every contig in the BAM. For each:
      1. Strip name_suffix to get the bare gene name.
      2. If restrict_genes is given and gene not in it -> skip.
      3. Determine utr5/utr3 (from table, or fixed values).
      4. Build the ±window profile anchored at true ORF start/end.
    Returns (tss_matrix, tes_matrix, n_skipped_restrict, n_skipped_na,
             n_skipped_short, n_total_contigs).
    """
    tss_rows, tes_rows = [], []
    n_total = n_skip_restrict = n_skip_na = n_skip_short = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = list(zip(bam.references, bam.lengths))
        n_total = len(refs)

        for i, (ref_name, ref_len) in enumerate(refs):
            if i % 1000 == 0:
                print(f"  {i}/{n_total}: {ref_name}")

            gene = strip_suffix(ref_name, name_suffix)

            if restrict_genes is not None and gene not in restrict_genes:
                n_skip_restrict += 1
                continue

            if utr_lookup is not None:
                utrs = utr_lookup.get(gene)
                if utrs is None:
                    n_skip_na += 1
                    continue
                utr5, utr3 = utrs
            else:
                utr5, utr3 = fixed_utr5, fixed_utr3

            orf_start = utr5
            orf_end   = ref_len - utr3
            orf_len   = orf_end - orf_start

            if orf_len < min_orf or orf_start < 0 or orf_end > ref_len:
                n_skip_short += 1
                continue

            try:
                cov = np.array(
                    bam.count_coverage(
                        ref_name,
                        quality_threshold=0,
                        read_callback=lambda r: (
                            not r.is_secondary and
                            not r.is_supplementary and
                            not r.is_unmapped
                        )
                    )
                ).sum(axis=0).astype(float)
            except (ValueError, KeyError):
                continue

            if len(cov) != ref_len:
                continue

            # TSS profile: position 0 = ORF start
            tss = np.zeros(2 * window, dtype=float)
            up_len = min(window, orf_start)
            if up_len > 0:
                tss[window - up_len: window] = cov[orf_start - up_len: orf_start]
            dn_len = min(window, orf_len)
            if dn_len > 0:
                tss[window: window + dn_len] = cov[orf_start: orf_start + dn_len]
            tss_rows.append(tss)

            # TES profile: position 0 = ORF end
            tes = np.zeros(2 * window, dtype=float)
            up_len = min(window, orf_len)
            if up_len > 0:
                tes[window - up_len: window] = cov[orf_end - up_len: orf_end]
            dn_len = min(window, ref_len - orf_end)
            if dn_len > 0:
                tes[window: window + dn_len] = cov[orf_end: orf_end + dn_len]
            tes_rows.append(tes)

    tss_mat = np.array(tss_rows) if tss_rows else np.zeros((0, 2 * window))
    tes_mat = np.array(tes_rows) if tes_rows else np.zeros((0, 2 * window))
    return tss_mat, tes_mat, n_skip_restrict, n_skip_na, n_skip_short, n_total


def write_profile(matrix, path, window):
    n = matrix.shape[0]
    if n == 0:
        print(f"WARNING: no transcripts -> {path} (file not written)")
        return
    mean = np.nanmean(matrix, axis=0)
    sem  = np.nanstd(matrix, axis=0) / np.sqrt(n)
    pos  = np.arange(-window, window)
    with open(path, "w") as fh:
        fh.write("position\tmean_coverage\tsem_coverage\tn_transcripts\n")
        for p, m, s in zip(pos, mean, sem):
            fh.write(f"{p}\t{m:.6f}\t{s:.6f}\t{n}\n")
    print(f"Written {n} transcripts -> {path}")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = Path(args.bam).stem
    print(f"\n=== {stem} ===")
    print(f"name_suffix='{args.name_suffix}'  window={args.window}  "
          f"min_orf={args.min_orf}")

    utr_lookup = None
    if args.utr_file:
        if not args.utr5_col or not args.utr3_col:
            sys.exit("--utr_file requires --utr5_col and --utr3_col")
        utr_lookup = load_utr_table(
            args.utr_file, args.utr5_col, args.utr3_col, args.gene_col)
        n_valid = sum(v is not None for v in utr_lookup.values())
        print(f"Loaded UTR table: {len(utr_lookup)} genes "
              f"({n_valid} with non-NA data)")
    else:
        print(f"Fixed UTRs: 5'={args.fixed_utr5}  3'={args.fixed_utr3}")

    restrict_genes = None
    if args.restrict_to:
        restrict_genes = load_restrict_set(args.restrict_to)
        print(f"Restricting to {len(restrict_genes)} genes "
              f"from {args.restrict_to}")

    tss_mat, tes_mat, n_skip_restrict, n_skip_na, n_skip_short, n_total = \
        compute_profiles(
            args.bam, utr_lookup, args.fixed_utr5, args.fixed_utr3,
            args.name_suffix, args.window, args.min_orf, restrict_genes)

    print(f"\nContigs in BAM:           {n_total}")
    print(f"Skipped (not in restrict): {n_skip_restrict}")
    print(f"Skipped (NA in UTR table): {n_skip_na}")
    print(f"Skipped (ORF too short / "
          f"bad coords):              {n_skip_short}")
    print(f"TSS transcripts used:      {tss_mat.shape[0]}")
    print(f"TES transcripts used:      {tes_mat.shape[0]}")

    if tss_mat.shape[0] == 0:
        print("\n*** WARNING: zero transcripts passed all filters. ***")
        print("Check that --name_suffix matches the BAM contig naming "
              "convention and that --restrict_to gene names are bare "
              "(no suffix) and match the UTR table's Gene column.")

    write_profile(tss_mat, outdir / f"{stem}_TSS.tsv", args.window)
    write_profile(tes_mat, outdir / f"{stem}_TES.tsv", args.window)
    print("Done.")


if __name__ == "__main__":
    main()
