#!/usr/bin/env python3
"""
derive_pelechano_utr.py
=======================
Derive per-gene 5'/3' UTR lengths from the Pelechano, Wei & Steinmetz (2013)
TIF-seq transcript-isoform table (GEO GSE39128, file GSE39128_tsedall.txt.gz),
so they can be compared against this study's final reference annotation in the
same way the Nagalakshmi et al. (2008) reference was compared.

Method
------
1. Read the ORF coordinates (the `gene` feature = ATG..stop genomic span,
   intron-inclusive) from the gene-models backbone GFF3. The two endpoints of
   the `gene` feature are exactly the outer CDS bounds, i.e. the ATG-side and
   stop-side genomic positions.
2. Read every TIF (transcript isoform; 5' end t5, 3' end t3, per-condition
   read counts) from the TIF-seq table.
3. For each ORF, select the *major* covering TIF (mTIF): among all TIFs on the
   same chromosome and strand whose t5 lies upstream of the ATG and whose t3
   lies downstream of the stop (i.e. that span the whole ORF), take the one
   with the highest total read support across all conditions.
4. Convert that mTIF's boundaries to UTR lengths:
       + strand : 5'UTR = ORF_start - t5 ; 3'UTR = t3 - ORF_end
       - strand : 5'UTR = t5 - ORF_end   ; 3'UTR = ORF_start - t3
   (UTR lengths are >= 0 by construction, since the mTIF covers the ORF.)
5. Genes with no covering TIF are written with NA UTR lengths, mirroring the
   NA convention of Nagalakshmi_UTR.csv.

Output
------
Pelechano_UTR.csv : CSV, columns Gene, five_prime_utr, three_prime_utr
                    (same schema as Nagalakshmi_UTR.csv)

Requirements: numpy, pandas
"""

import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path

# Roman-numeral chromosome names (GFF3) -> Arabic (TIF-seq table uses 1..16)
ROMAN_TO_ARABIC = {
    "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, "VIII": 8,
    "IX": 9, "X": 10, "XI": 11, "XII": 12, "XIII": 13, "XIV": 14, "XV": 15,
    "XVI": 16,
}

COUNT_COLS = ["ypd", "gal", "lypd", "lgal", "nypd", "ngal"]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff3", required=True,
                   help="gene-models backbone GFF3 (provides ORF/gene spans)")
    p.add_argument("--tif", required=True,
                   help="TIF-seq table (GSE39128_tsedall.txt.gz, gz or plain)")
    p.add_argument("--out", required=True,
                   help="output Pelechano_UTR.csv path")
    return p.parse_args()


def load_orfs(gff3_path):
    """Read `gene` features -> DataFrame(gene, chr, strand, start, end)."""
    rows = []
    id_re = re.compile(r"ID=([^;]+)")
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = ROMAN_TO_ARABIC.get(f[0].strip())
            if chrom is None:
                continue  # skip mito / unrecognised contigs
            try:
                start = int(f[3]); end = int(f[4])
            except ValueError:
                continue
            m = id_re.search(f[8])
            if not m:
                continue
            rows.append((m.group(1), chrom, f[6], start, end))
    orfs = pd.DataFrame(rows, columns=["Gene", "chr", "strand", "start", "end"])
    print(f"ORFs loaded: {len(orfs)} (chr 1-16, gene features)")
    return orfs


def load_tifs(tif_path):
    """Read TIF table; add total read support across conditions."""
    df = pd.read_csv(tif_path, sep="\t")
    df["total"] = df[COUNT_COLS].sum(axis=1)
    # keep only the columns we need
    df = df[["chr", "strand", "t5", "t3", "total"]].copy()
    print(f"TIFs loaded: {len(df):,} rows")
    return df


def derive(orfs, tifs):
    """For each ORF, pick the highest-support covering TIF and compute UTRs."""
    out_five = np.full(len(orfs), np.nan)
    out_three = np.full(len(orfs), np.nan)
    out_support = np.zeros(len(orfs), dtype=int)

    # Pre-split TIFs by (chr, strand) into numpy arrays for fast filtering.
    groups = {}
    for (c, s), g in tifs.groupby(["chr", "strand"]):
        groups[(c, s)] = (g["t5"].to_numpy(), g["t3"].to_numpy(),
                          g["total"].to_numpy())

    n_covered = 0
    for i, row in enumerate(orfs.itertuples(index=False)):
        key = (row.chr, row.strand)
        if key not in groups:
            continue
        t5, t3, tot = groups[key]
        if row.strand == "+":
            # transcript spans ORF: t5 upstream of ATG, t3 downstream of stop
            mask = (t5 <= row.start) & (t3 >= row.end)
        else:
            # minus strand: 5' end = high coord, 3' end = low coord
            mask = (t5 >= row.end) & (t3 <= row.start)
        if not mask.any():
            continue
        idx = np.flatnonzero(mask)
        best = idx[np.argmax(tot[idx])]
        b5, b3 = t5[best], t3[best]
        if row.strand == "+":
            five = row.start - b5
            three = b3 - row.end
        else:
            five = b5 - row.end
            three = row.start - b3
        out_five[i] = max(0, five)
        out_three[i] = max(0, three)
        out_support[i] = tot[best]
        n_covered += 1

    res = pd.DataFrame({
        "Gene": orfs["Gene"].values,
        "five_prime_utr": out_five,
        "three_prime_utr": out_three,
        "mtif_support": out_support,
    })
    print(f"Genes with a covering mTIF: {n_covered} / {len(orfs)}")
    return res


def main():
    args = parse_args()
    orfs = load_orfs(args.gff3)
    tifs = load_tifs(args.tif)
    res = derive(orfs, tifs)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Write the Nagalakshmi-style 3-column table (Gene, 5'UTR, 3'UTR);
    # keep mtif_support in a sidecar for QC.
    res[["Gene", "five_prime_utr", "three_prime_utr"]].to_csv(out, index=False)
    res.to_csv(out.with_suffix(".qc.csv"), index=False)
    print(f"Written: {out}")
    print(f"  median 5'UTR = {res['five_prime_utr'].median():.0f} nt, "
          f"median 3'UTR = {res['three_prime_utr'].median():.0f} nt "
          f"(over genes with a covering mTIF)")


if __name__ == "__main__":
    main()
