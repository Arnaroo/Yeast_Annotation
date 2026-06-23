#!/usr/bin/env python3
"""
plot_read_pileup_multigene_with_nag_drs.py
---------------------------------------------
Same as plot_read_pileup_multigene.py (read pileup, one panel per gene,
single BAM), but additionally overlays:
  - Nagalakshmi's annotation as a second (purple) track, drawn below
    our annotation track
  - DRS-derived (pre-merge) UTR boundaries as dashed vertical lines in
    the read panel, in a third colour, distinct from our solid black
    ORF-boundary lines and Nagalakshmi's own boundary markers

This is a separate script — plot_read_pileup_multigene.py is untouched.

Coordinate system: x = 0 is always OUR annotated ORF start (same as
plot_read_pileup_multigene.py). Nagalakshmi and DRS boundaries are
placed in this same coordinate system:
  Nagalakshmi 5' boundary: x = -nag_utr5
  Nagalakshmi 3' boundary: x = orf_len + nag_utr3
  DRS 5' boundary:         x = -drs_utr5
  DRS 3' boundary:         x = orf_len + drs_utr3
(orf_len is computed from OUR annotation; the ORF itself is identical
across all three annotations since they all annotate the same SGD ORF.)

Usage:
    python plot_read_pileup_multigene_with_nag_drs.py \\
        --bam          sorted_reads.bam \\
        --genes        GENE_A GENE_B \\
        --panel_labels A B \\
        --utr_our      final_utr.tsv \\
        --utr_nag      reference_UTR.csv \\
        --utr_drs      DRS_UTR_corrected.tsv \\
        --name_suffix  _mRNA \\
        --outdir       ./figures \\
        --outname      pileup_drs_vs_ref

Requirements: pysam, numpy, matplotlib, pandas
"""

import argparse
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

C_FWD = "#4477AA"
C_REV = "#EE6677"
C_OUR = "#228833"
C_NAG = "#AA3377"
C_DRS = "#CC6600"

FS_BODY  = 8
FS_LABEL = 9
FS_PANEL = 12

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FS_BODY,
    "axes.labelsize":    FS_LABEL,
    "xtick.labelsize":   FS_BODY,
    "ytick.labelsize":   FS_BODY,
    "axes.linewidth":    0.8,
    "xtick.major.width": 0.8,
})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam",          required=True,
                   help="Sorted, indexed BAM file")
    p.add_argument("--genes",        nargs="+", required=True,
                   help="Gene names (must match BAM contig names + name_suffix)")
    p.add_argument("--panel_labels", nargs="+", default=None,
                   help="Panel letters (default: A, B, C, …)")
    p.add_argument("--utr_our",      required=True,
                   help="Final merged UTR annotation TSV "
                        "(columns: Gene  final_five_prime_utr  final_three_prime_utr)")
    p.add_argument("--utr_nag",      required=True,
                   help="Reference UTR annotation CSV "
                        "(columns: Gene  five_prime_utr  three_prime_utr; NAs allowed)")
    p.add_argument("--utr_drs",      required=True,
                   help="Pre-merge DRS UTR TSV "
                        "(columns: gene  max_five_prime_utr  max_three_prime_utr). "
                        "Must be the genuine pre-merge file, not a copy of the "
                        "final merged annotation.")
    p.add_argument("--name_suffix",  default="_mRNA",
                   help="Suffix appended to gene name to form BAM contig (default: _mRNA)")
    p.add_argument("--outdir",       required=True,
                   help="Output directory")
    p.add_argument("--outname",      default="pileup_drs_vs_ref",
                   help="Output filename stem")
    p.add_argument("--max_reads",    type=int, default=300)
    p.add_argument("--seed",         type=int, default=0)
    return p.parse_args()


def load_utr_our(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = ["Gene", "utr5", "utr3"]
    return df.set_index("Gene")


def load_utr_nag(path):
    df = pd.read_csv(path)
    df.columns = ["Gene", "utr5", "utr3"]
    df["utr5"] = pd.to_numeric(df["utr5"], errors="coerce")
    df["utr3"] = pd.to_numeric(df["utr3"], errors="coerce")
    return df.set_index("Gene")


def load_utr_drs(path):
    """
    Load the pre-merge DRS-derived UTR table (max_utr.tsv / 
    DRS_UTR_corrected.tsv), NOT the mislabeled DRS_UTR.tsv which is
    actually a copy of the post-merge final_utr.tsv.
    Expected columns: gene, max_five_prime_utr, max_three_prime_utr
    (3 columns, no cds column — this table only has ~5,417 genes,
    fewer than our_ref/Nagalakshmi, since genes without DRS expression
    are simply absent).
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = ["Gene", "utr5", "utr3"]
    return df.set_index("Gene")


def get_ref_len(bam_path, contig):
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        d = dict(zip(bam.references, bam.lengths))
        return d.get(contig)


def aligned_blocks(read):
    """Aligned (start, end) blocks on the reference, split at N gaps."""
    blocks = []
    ref_pos = read.reference_start
    cur_start = ref_pos
    cur_len = 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            cur_len += length
            ref_pos += length
        elif op == 2:
            cur_len += length
            ref_pos += length
        elif op == 3:
            blocks.append((cur_start, cur_start + cur_len))
            ref_pos += length
            cur_start = ref_pos
            cur_len = 0
        elif op in (1, 4, 5, 6):
            continue
    if cur_len > 0:
        blocks.append((cur_start, cur_start + cur_len))
    return blocks


def fetch_reads(bam_path, contig, ref_len, max_reads, seed):
    reads = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        try:
            for r in bam.fetch(contig, 0, ref_len):
                if r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                blocks = aligned_blocks(r)
                if not blocks:
                    continue
                reads.append({
                    "blocks": blocks,
                    "strand": "-" if r.is_reverse else "+",
                    "start":  blocks[0][0],
                    "end":    blocks[-1][1],
                })
        except ValueError:
            return reads
    if len(reads) > max_reads:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(reads), size=max_reads, replace=False)
        reads = [reads[i] for i in sorted(idx)]
    return reads


def assign_rows(reads, gap=15):
    reads_sorted = sorted(reads, key=lambda r: r["start"])
    row_ends = []
    for r in reads_sorted:
        placed = False
        for i, end in enumerate(row_ends):
            if r["start"] > end + gap:
                r["row"] = i
                row_ends[i] = r["end"]
                placed = True
                break
        if not placed:
            r["row"] = len(row_ends)
            row_ends.append(r["end"])
    return reads_sorted, len(row_ends)


def draw_reads(ax, reads, n_rows, utr5, height=0.8):
    for r in reads:
        color = C_FWD if r["strand"] == "+" else C_REV
        y = r["row"]
        x0 = r["blocks"][0][0] - utr5
        x1 = r["blocks"][-1][1] - utr5
        ax.plot([x0, x1], [y + height / 2] * 2, color=color, lw=0.6,
                alpha=0.5, zorder=1)
        for b_start, b_end in r["blocks"]:
            ax.add_patch(mpatches.Rectangle(
                (b_start - utr5, y), b_end - b_start, height,
                facecolor=color, edgecolor="none", alpha=0.75, zorder=2))
    ax.set_ylim(-0.5, max(n_rows, 1) + 0.5)
    ax.invert_yaxis()
    ax.set_yticks([])


def draw_annotation_dual(ax, utr5_our, utr3_our, orf_len,
                         utr5_nag, utr3_nag, gene, height=0.62):
    """
    Two annotation rows in one strip:
      y=0.00: our annotation (green) — gene name embedded as white text
      y=1.10: Nagalakshmi (purple)   — label embedded as white text
    """
    def draw_track(y0, u5, u3, color, label):
        if u5 and u5 > 0:
            ax.add_patch(mpatches.Rectangle(
                (-u5, y0), u5, height,
                facecolor=color, alpha=0.30, edgecolor=color, lw=0.8))
        ax.add_patch(mpatches.Rectangle(
            (0, y0), orf_len, height,
            facecolor=color, alpha=0.85, edgecolor=color, lw=0.8))
        if u3 and u3 > 0:
            ax.add_patch(mpatches.Rectangle(
                (orf_len, y0), u3, height,
                facecolor=color, alpha=0.30, edgecolor=color, lw=0.8))
        # Label centred in the ORF block, white bold text
        ax.text(orf_len / 2, y0 + height / 2, label,
                ha="center", va="center",
                fontsize=FS_BODY, fontweight="bold",
                color="white", zorder=10, clip_on=True)

    draw_track(0.0, utr5_our, utr3_our, C_OUR,
               f"Our ref  ·  {gene}")
    if utr5_nag is not None:
        draw_track(1.10, utr5_nag, utr3_nag, C_NAG,
                   "Nagalakshmi")


def prepare_gene_data(gene, bam_path, utr_our_df, utr_nag_df, utr_drs_df,
                      name_suffix, max_reads, seed):
    contig = gene + name_suffix

    if gene not in utr_our_df.index:
        return {"gene": gene, "error": f"{gene}: not in our UTR table"}

    utr5_our = int(utr_our_df.loc[gene, "utr5"])
    utr3_our = int(utr_our_df.loc[gene, "utr3"])

    utr5_nag = utr3_nag = None
    if gene in utr_nag_df.index:
        u5n = utr_nag_df.loc[gene, "utr5"]
        u3n = utr_nag_df.loc[gene, "utr3"]
        if pd.notna(u5n) and pd.notna(u3n):
            utr5_nag, utr3_nag = int(u5n), int(u3n)

    utr5_drs = utr3_drs = None
    if gene in utr_drs_df.index:
        u5d = utr_drs_df.loc[gene, "utr5"]
        u3d = utr_drs_df.loc[gene, "utr3"]
        if pd.notna(u5d) and pd.notna(u3d):
            utr5_drs, utr3_drs = int(u5d), int(u3d)

    ref_len = get_ref_len(bam_path, contig)
    if ref_len is None:
        return {"gene": gene,
                "error": f"{gene}: contig '{contig}' not in BAM"}

    orf_len = ref_len - utr5_our - utr3_our
    reads = fetch_reads(bam_path, contig, ref_len, max_reads, seed)
    reads, n_rows = assign_rows(reads)
    print(f"{gene}: our 5'UTR={utr5_our} ORF={orf_len} 3'UTR={utr3_our} | "
          f"Nag 5'={utr5_nag} 3'={utr3_nag} | "
          f"DRS 5'={utr5_drs} 3'={utr3_drs} | "
          f"{len(reads)} reads, {n_rows} rows")

    return {
        "gene": gene, "error": None, "contig": contig,
        "utr5_our": utr5_our, "utr3_our": utr3_our, "orf_len": orf_len,
        "utr5_nag": utr5_nag, "utr3_nag": utr3_nag,
        "utr5_drs": utr5_drs, "utr3_drs": utr3_drs,
        "reads": reads, "n_rows": n_rows,
    }


def draw_panel(ax_annot, ax_reads, panel_label, data):
    gene = data["gene"]

    if data["error"]:
        ax_reads.text(0.5, 0.5, data["error"], ha="center", va="center",
                      transform=ax_reads.transAxes, color="grey")
        ax_annot.set_yticks([]); ax_annot.set_xticks([])
        for spine in ax_annot.spines.values():
            spine.set_visible(False)
        ax_annot.text(-0.01, 0.5, panel_label, transform=ax_annot.transAxes,
                      fontsize=13, fontweight="bold", ha="right", va="center")
        return

    utr5_our, utr3_our = data["utr5_our"], data["utr3_our"]
    orf_len = data["orf_len"]
    utr5_nag, utr3_nag = data["utr5_nag"], data["utr3_nag"]
    utr5_drs, utr3_drs = data["utr5_drs"], data["utr3_drs"]
    reads, n_rows = data["reads"], data["n_rows"]

    # Dual annotation strip (our + Nagalakshmi); gene name inside Our ref bar
    draw_annotation_dual(ax_annot, utr5_our, utr3_our, orf_len,
                         utr5_nag, utr3_nag, gene)
    ax_annot.set_ylim(1.9, -0.3)
    ax_annot.set_yticks([]); ax_annot.set_xticks([])
    for spine in ax_annot.spines.values():
        spine.set_visible(False)

    ax_annot.text(-0.01, 0.5, panel_label, transform=ax_annot.transAxes,
                  fontsize=FS_PANEL, fontweight="bold", ha="right", va="center")

    # Reads
    draw_reads(ax_reads, reads, n_rows, utr5_our)
    ax_reads.axvspan(0, orf_len, alpha=0.05, color="steelblue", zorder=0)

    # Our ORF boundaries — solid black dashed
    ax_reads.axvline(0,       color="black", lw=0.9, ls="--", alpha=0.55,
                     label="_nolegend_")
    ax_reads.axvline(orf_len, color="black", lw=0.9, ls="--", alpha=0.55,
                     label="_nolegend_")

    # Nagalakshmi boundaries — purple dotted
    if utr5_nag is not None:
        ax_reads.axvline(-utr5_nag, color=C_NAG, lw=1.0, ls=":", alpha=0.8,
                         label="_nolegend_")
        ax_reads.axvline(orf_len + utr3_nag, color=C_NAG, lw=1.0, ls=":",
                         alpha=0.8, label="_nolegend_")

    # DRS boundaries — orange dashed (distinct from both above)
    if utr5_drs is not None:
        ax_reads.axvline(-utr5_drs, color=C_DRS, lw=1.0, ls=(0, (4, 2)),
                         alpha=0.85, label="_nolegend_")
        ax_reads.axvline(orf_len + utr3_drs, color=C_DRS, lw=1.0,
                         ls=(0, (4, 2)), alpha=0.85, label="_nolegend_")

    ax_reads.spines["top"].set_visible(False)
    ax_reads.spines["right"].set_visible(False)
    ax_reads.spines["left"].set_visible(False)

    info = (
        f"Our: 5\u2032={utr5_our}  3\u2032={utr3_our}"
        f"   \u2502   Nag: 5\u2032={utr5_nag}  3\u2032={utr3_nag}"
        f"   \u2502   DRS (pre-merge): 5\u2032={utr5_drs}  3\u2032={utr3_drs}"
        f"   \u2502   n = {len(reads)} reads"
    )
    ax_reads.text(0.01, 0.03, info, transform=ax_reads.transAxes,
                  fontsize=7.5, fontweight="bold", va="bottom", ha="left",
                  bbox=dict(boxstyle="round,pad=0.2", fc="white",
                            ec="none", alpha=0.75))

    xmin = -max(utr5_our, utr5_nag or 0, utr5_drs or 0) - 100
    xmax = orf_len + max(utr3_our, utr3_nag or 0, utr3_drs or 0) + 100
    ax_annot.set_xlim(xmin, xmax)
    ax_reads.set_xlim(xmin, xmax)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = args.genes
    labels = args.panel_labels or [chr(ord("A") + i) for i in range(len(genes))]
    if len(labels) != len(genes):
        raise SystemExit("--panel_labels must match --genes in count")

    utr_our_df = load_utr_our(args.utr_our)
    utr_nag_df = load_utr_nag(args.utr_nag)
    utr_drs_df = load_utr_drs(args.utr_drs)

    gene_data = [
        prepare_gene_data(g, args.bam, utr_our_df, utr_nag_df, utr_drs_df,
                          args.name_suffix, args.max_reads, args.seed)
        for g in genes
    ]

    ROW_HEIGHT_IN = 0.045
    MIN_READS_H   = 1.8
    MAX_READS_H   = 7.0
    ANNOT_H       = 0.65     # taller — two rows with embedded text
    GENE_GAP      = 0.55

    height_ratios, panel_kind = [], []
    for i, d in enumerate(gene_data):
        n_rows = d.get("n_rows", 3) if not d["error"] else 3
        reads_h = min(max(n_rows * ROW_HEIGHT_IN, MIN_READS_H), MAX_READS_H)
        height_ratios += [ANNOT_H, reads_h]
        panel_kind += ["annot", "reads"]
        if i < len(gene_data) - 1:
            height_ratios.append(GENE_GAP)
            panel_kind.append("gap")

    fig_height = sum(height_ratios) + 0.7
    fig = plt.figure(figsize=(9, fig_height))
    gs = fig.add_gridspec(len(height_ratios), 1,
                          height_ratios=height_ratios, hspace=0)

    axes_by_gene = []
    row_i = 0
    for kind in panel_kind:
        ax = fig.add_subplot(gs[row_i])
        if kind == "gap":
            ax.axis("off")
        else:
            axes_by_gene.append(ax)
        row_i += 1

    for i, (gene, label, d) in enumerate(zip(genes, labels, gene_data)):
        ax_annot = axes_by_gene[i * 2]
        ax_reads = axes_by_gene[i * 2 + 1]
        draw_panel(ax_annot, ax_reads, label, d)
        ax_reads.set_xlabel("Position relative to ORF start (bp)", fontsize=9)

    fwd_patch = mpatches.Patch(color=C_FWD, label="Forward strand read")
    rev_patch = mpatches.Patch(color=C_REV, label="Reverse strand read")
    our_patch = mpatches.Patch(color=C_OUR, label="Our annotation")
    nag_patch = mpatches.Patch(color=C_NAG, label="Nagalakshmi annotation")
    drs_line  = plt.Line2D([0], [0], color=C_DRS, lw=1.5, ls=(0, (4, 2)),
                           label="DRS boundary (pre-merge)")
    fig.legend(handles=[fwd_patch, rev_patch, our_patch, nag_patch, drs_line],
               loc="upper right", fontsize=FS_BODY, frameon=False,
               bbox_to_anchor=(0.99, 0.998), ncol=1,
               handlelength=1.6, handletextpad=0.5)

    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.subplots_adjust(hspace=0)

    for ext in ("svg", "png", "pdf"):
        outpath = outdir / f"{args.outname}.{ext}"
        dpi = 300 if ext == "png" else None
        plt.savefig(outpath, bbox_inches="tight", dpi=dpi)
        print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
