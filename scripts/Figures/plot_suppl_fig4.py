#!/usr/bin/env python3
"""
plot_suppl_fig4.py
-------------------
Supplementary Figure 4: read pileup for DRS and short-read (SR) BAMs,
one gene per panel group, with our annotation + Nagalakshmi annotation
tracks and pre-merge DRS boundary lines.

Layout per gene:
  ┌─────────────────────────────────────────┐
  │ annotation strip (Our ref + Nagalakshmi)│
  ├─────────────────────────────────────────┤
  │ DRS read pileup                         │
  ├─────────────────────────────────────────┤
  │ SR read pileup                          │
  └─────────────────────────────────────────┘

All panels share the same x-axis (position relative to OUR ORF start).
Vertical lines in read panels:
  black dashed  : our ORF boundaries (start=0, end=orf_len)
  purple dotted : Nagalakshmi UTR boundaries (-nag_utr5, orf_len+nag_utr3)
  orange dashed : pre-merge DRS boundaries (-drs_utr5, orf_len+drs_utr3)
Info label bottom-left of each read panel (not obscuring reads).

Usage:
    python plot_suppl_fig4.py \\
        --drs_bam      DRS_reads.sorted.bam \\
        --sr_bam       SR_reads.sorted.bam \\
        --genes        GENE_A GENE_B GENE_C GENE_D \\
        --panel_labels A B C D \\
        --bam1_label   "DRS (independent)" \\
        --bam2_label   "Short Read (independent)" \\
        --utr_our      final_utr.tsv \\
        --utr_nag      reference_UTR.csv \\
        --utr_drs      DRS_UTR_corrected.tsv \\
        --name_suffix  _mRNA \\
        --outdir       ./figures \\
        --outname      suppl_fig4

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

# ── Colours ───────────────────────────────────────────────────────────────
C_DRS1_FWD = "#4477AA"   # DRS forward reads (blue)
C_DRS1_REV = "#88CCEE"   # DRS reverse reads (light blue)
C_DRS2_FWD = "#EE6677"   # Short-read forward reads (salmon)
C_DRS2_REV = "#FFAABB"   # Short-read reverse reads (light salmon)
C_OUR      = "#228833"   # Our annotation (green)
C_NAG      = "#AA3377"   # Nagalakshmi annotation (purple)
C_DRS_BND  = "#CC6600"   # Pre-merge DRS boundary lines (orange)

# Aliases kept for compatibility with draw_reads calls
C_DRS_READ = C_DRS1_FWD
C_DRS_REV  = C_DRS1_REV
C_SR_READ  = C_DRS2_FWD
C_SR_REV   = C_DRS2_REV

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
    p.add_argument("--drs_bam",      required=True,
                   help="Sorted, indexed DRS (long-read) BAM file, "
                        "mapped to the same reference as --utr_our")
    p.add_argument("--sr_bam",       required=True,
                   help="Sorted, indexed short-read BAM file, "
                        "mapped to the same reference as --utr_our")
    p.add_argument("--genes",        nargs="+", required=True,
                   help="Exactly 4 gene names (2×2 panel grid)")
    p.add_argument("--panel_labels", nargs="+", default=None,
                   help="Panel letters (default: A, B, C, D)")
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
    p.add_argument("--bam1_label",   default="DRS (independent)",
                   help="Legend label for the DRS BAM panel")
    p.add_argument("--bam2_label",   default="Short Read (independent)",
                   help="Legend label for the short-read BAM panel")
    p.add_argument("--outdir",       required=True,
                   help="Output directory")
    p.add_argument("--outname",      default="suppl_fig4",
                   help="Output filename stem")
    p.add_argument("--max_reads",    type=int, default=300)
    p.add_argument("--seed",         type=int, default=0)
    return p.parse_args()


# ── UTR loaders ───────────────────────────────────────────────────────────

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
    Pre-merge DRS UTR table (DRS_UTR_corrected.tsv / max_utr.tsv).
    Columns: gene, max_five_prime_utr, max_three_prime_utr  (3 cols).
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = ["Gene", "utr5", "utr3"]
    return df.set_index("Gene")


# ── Read helpers ──────────────────────────────────────────────────────────

def get_ref_len(bam_path, contig):
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return dict(zip(bam.references, bam.lengths)).get(contig)


def aligned_blocks(read):
    blocks, ref_pos, cur_start, cur_len = [], read.reference_start, read.reference_start, 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            cur_len += length; ref_pos += length
        elif op == 2:
            cur_len += length; ref_pos += length
        elif op == 3:
            blocks.append((cur_start, cur_start + cur_len))
            ref_pos += length; cur_start = ref_pos; cur_len = 0
        # I, S, H, P: skip
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
                reads.append({"blocks": blocks,
                              "strand": "-" if r.is_reverse else "+",
                              "start": blocks[0][0],
                              "end":   blocks[-1][1]})
        except ValueError:
            pass
    if len(reads) > max_reads:
        idx = np.random.default_rng(seed).choice(len(reads), max_reads,
                                                  replace=False)
        reads = [reads[i] for i in sorted(idx)]
    return reads


def assign_rows(reads, gap=15):
    reads = sorted(reads, key=lambda r: r["start"])
    row_ends = []
    for r in reads:
        placed = False
        for i, end in enumerate(row_ends):
            if r["start"] > end + gap:
                r["row"] = i; row_ends[i] = r["end"]; placed = True; break
        if not placed:
            r["row"] = len(row_ends); row_ends.append(r["end"])
    return reads, len(row_ends)


# ── Drawing helpers ───────────────────────────────────────────────────────

def draw_reads(ax, reads, n_rows, utr5, c_fwd, c_rev, height=0.8):
    for r in reads:
        color = c_fwd if r["strand"] == "+" else c_rev
        y = r["row"]
        x0 = r["blocks"][0][0] - utr5
        x1 = r["blocks"][-1][1] - utr5
        ax.plot([x0, x1], [y + height / 2] * 2,
                color=color, lw=0.5, alpha=0.4, zorder=1)
        for b_start, b_end in r["blocks"]:
            ax.add_patch(mpatches.Rectangle(
                (b_start - utr5, y), b_end - b_start, height,
                facecolor=color, edgecolor="none", alpha=0.75, zorder=2))
    ax.set_ylim(-0.5, max(n_rows, 1) + 0.5)
    ax.invert_yaxis()
    ax.set_yticks([])


def draw_vlines(ax, orf_len, utr5_nag, utr3_nag, utr5_drs, utr3_drs):
    """Draw ORF, Nagalakshmi and DRS boundary lines on a read panel."""
    ax.axvspan(0, orf_len, alpha=0.05, color="steelblue", zorder=0)
    ax.axvline(0,       color="black", lw=0.9, ls="--", alpha=0.5)
    ax.axvline(orf_len, color="black", lw=0.9, ls="--", alpha=0.5)
    if utr5_nag is not None:
        ax.axvline(-utr5_nag,         color=C_NAG, lw=1.0, ls=":",
                   alpha=0.8)
        ax.axvline(orf_len + utr3_nag, color=C_NAG, lw=1.0, ls=":",
                   alpha=0.8)
    if utr5_drs is not None:
        ax.axvline(-utr5_drs,         color=C_DRS_BND, lw=1.0,
                   ls=(0, (4, 2)), alpha=0.85)
        ax.axvline(orf_len + utr3_drs, color=C_DRS_BND, lw=1.0,
                   ls=(0, (4, 2)), alpha=0.85)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def draw_info(ax, label, utr5_our, utr3_our, orf_len,
              utr5_nag, utr3_nag, utr5_drs, utr3_drs, n_reads):
    nag_str = (f"Nag: 5\u2032={utr5_nag}  3\u2032={utr3_nag}"
               if utr5_nag is not None else "Nag: \u2014")
    drs_str = (f"DRS: 5\u2032={utr5_drs}  3\u2032={utr3_drs}"
               if utr5_drs is not None else "DRS: \u2014")
    text = (f"Our: 5\u2032={utr5_our}  3\u2032={utr3_our}"
            f"   \u2502   {nag_str}"
            f"   \u2502   {drs_str}"
            f"   \u2502   {label} (n={n_reads})")
    ax.text(0.01, 0.03, text, transform=ax.transAxes, fontsize=6.5,
            fontweight="bold", va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


def draw_annot_strip(ax, utr5_our, utr3_our, orf_len,
                     utr5_nag, utr3_nag, gene, height=0.62):
    """Two annotation rows; gene name and track labels embedded as white text."""
    def track(y0, u5, u3, color, label):
        if u5 and u5 > 0:
            ax.add_patch(mpatches.Rectangle((-u5, y0), u5, height,
                         facecolor=color, alpha=0.30, edgecolor=color, lw=0.8))
        ax.add_patch(mpatches.Rectangle((0, y0), orf_len, height,
                     facecolor=color, alpha=0.85, edgecolor=color, lw=0.8))
        if u3 and u3 > 0:
            ax.add_patch(mpatches.Rectangle((orf_len, y0), u3, height,
                         facecolor=color, alpha=0.30, edgecolor=color, lw=0.8))
        ax.text(orf_len / 2, y0 + height / 2, label,
                ha="center", va="center", fontsize=FS_BODY,
                fontweight="bold", color="white", zorder=10, clip_on=True)

    track(0.0, utr5_our, utr3_our, C_OUR, f"Our ref  \u00b7  {gene}")
    if utr5_nag is not None:
        track(1.1, utr5_nag, utr3_nag, C_NAG, "Nagalakshmi")
    ax.set_ylim(1.9, -0.3)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


# ── Per-gene data preparation ─────────────────────────────────────────────

def prepare_gene(gene, drs_bam, sr_bam, utr_our_df, utr_nag_df,
                 utr_drs_df, suffix, max_reads, seed):
    contig = gene + suffix
    error = None

    if gene not in utr_our_df.index:
        return {"gene": gene, "error": f"{gene}: not in our UTR table"}

    utr5_our = int(utr_our_df.loc[gene, "utr5"])
    utr3_our = int(utr_our_df.loc[gene, "utr3"])

    def _nag(col):
        if gene in utr_nag_df.index:
            v = utr_nag_df.loc[gene, col]
            return int(v) if pd.notna(v) else None
        return None

    utr5_nag = _nag("utr5"); utr3_nag = _nag("utr3")

    def _drs(col):
        if gene in utr_drs_df.index:
            v = utr_drs_df.loc[gene, col]
            return int(v) if pd.notna(v) else None
        return None

    utr5_drs = _drs("utr5"); utr3_drs = _drs("utr3")

    ref_len = get_ref_len(drs_bam, contig)
    if ref_len is None:
        ref_len = get_ref_len(sr_bam, contig)
    if ref_len is None:
        return {"gene": gene,
                "error": f"{gene}: contig '{contig}' not in either BAM"}

    orf_len = ref_len - utr5_our - utr3_our

    drs_reads, drs_rows = assign_rows(
        fetch_reads(drs_bam, contig, ref_len, max_reads, seed))
    sr_reads,  sr_rows  = assign_rows(
        fetch_reads(sr_bam,  contig, ref_len, max_reads, seed))

    print(f"{gene}: 5'={utr5_our} ORF={orf_len} 3'={utr3_our} | "
          f"Nag 5'={utr5_nag} 3'={utr3_nag} | DRS 5'={utr5_drs} 3'={utr3_drs} | "
          f"DRS reads={len(drs_reads)} rows={drs_rows} | "
          f"SR reads={len(sr_reads)} rows={sr_rows}")

    return {
        "gene": gene, "error": None, "contig": contig,
        "utr5_our": utr5_our, "utr3_our": utr3_our, "orf_len": orf_len,
        "utr5_nag": utr5_nag, "utr3_nag": utr3_nag,
        "utr5_drs": utr5_drs, "utr3_drs": utr3_drs,
        "drs_reads": drs_reads, "drs_rows": drs_rows,
        "sr_reads":  sr_reads,  "sr_rows":  sr_rows,
    }


# ── Figure builder (2 × 2 layout) ────────────────────────────────────────

def make_figure(genes, labels, gene_data, outdir, outname,
                bam1_label="DRS (NS condition)",
                bam2_label="DRS (S10 condition)"):
    """
    2-row × 2-col layout: each cell contains one gene with
    annotation strip / DRS pileup / SR pileup stacked vertically.
    Row heights are equalised within each row for clean alignment.
    """
    assert len(genes) == 4, "Exactly 4 genes required for 2×2 layout"

    ROW_H   = 0.045   # inches per stacked read row
    MIN_H   = 1.8
    MAX_H   = 4.5     # cap tighter than standalone to fit 2×2
    ANNOT_H = 0.65
    FIG_W   = 14.0

    # Compute DRS and SR panel heights per main-grid row
    # (each row contains 2 genes; heights equalised to taller one)
    def row_heights(idx_a, idx_b):
        dh = max(
            min(max(gene_data[idx_a].get("drs_rows", 3) * ROW_H, MIN_H), MAX_H),
            min(max(gene_data[idx_b].get("drs_rows", 3) * ROW_H, MIN_H), MAX_H),
        )
        sh = max(
            min(max(gene_data[idx_a].get("sr_rows", 3) * ROW_H, MIN_H), MAX_H),
            min(max(gene_data[idx_b].get("sr_rows", 3) * ROW_H, MIN_H), MAX_H),
        )
        return dh, sh

    drs_h0, sr_h0 = row_heights(0, 1)
    drs_h1, sr_h1 = row_heights(2, 3)

    row0_h = ANNOT_H + drs_h0 + sr_h0
    row1_h = ANNOT_H + drs_h1 + sr_h1
    fig_h  = row0_h + row1_h + 1.2   # 1.2 in for legend + margins

    fig = plt.figure(figsize=(FIG_W, fig_h), facecolor="white")

    # Main 2×2 gridspec
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[row0_h, row1_h],
        hspace=0.28, wspace=0.22,
        left=0.06, right=0.97,
        top=0.97, bottom=0.08,
    )

    all_drs_h = [drs_h0, drs_h0, drs_h1, drs_h1]
    all_sr_h  = [sr_h0,  sr_h0,  sr_h1,  sr_h1]

    for gi, (gene, label, d) in enumerate(zip(genes, labels, gene_data)):
        row, col = divmod(gi, 2)
        drs_h = all_drs_h[gi]
        sr_h  = all_sr_h[gi]

        # Nested subgridspec: annot / DRS / SR
        gs_cell = gs[row, col].subgridspec(
            3, 1,
            height_ratios=[ANNOT_H, drs_h, sr_h],
            hspace=0.0,
        )
        ax_annot = fig.add_subplot(gs_cell[0])
        ax_drs   = fig.add_subplot(gs_cell[1])
        ax_sr    = fig.add_subplot(gs_cell[2])

        if d["error"]:
            for ax in (ax_annot, ax_drs, ax_sr):
                ax.text(0.5, 0.5, d["error"], ha="center", va="center",
                        transform=ax.transAxes, color="grey", fontsize=FS_BODY)
                ax.axis("off")
            continue

        u5o, u3o = d["utr5_our"], d["utr3_our"]
        orf_len  = d["orf_len"]
        u5n, u3n = d["utr5_nag"], d["utr3_nag"]
        u5d, u3d = d["utr5_drs"], d["utr3_drs"]

        xmin = -(max(u5o, u5n or 0, u5d or 0) + 100)
        xmax = orf_len + max(u3o, u3n or 0, u3d or 0) + 100

        # Annotation strip
        draw_annot_strip(ax_annot, u5o, u3o, orf_len, u5n, u3n, gene)
        ax_annot.set_xlim(xmin, xmax)
        ax_annot.text(-0.03, 0.5, label, transform=ax_annot.transAxes,
                      fontsize=FS_PANEL, fontweight="bold",
                      ha="right", va="center")

        # DRS reads
        draw_reads(ax_drs, d["drs_reads"], d["drs_rows"], u5o,
                   C_DRS1_FWD, C_DRS1_REV)
        draw_vlines(ax_drs, orf_len, u5n, u3n, u5d, u3d)
        draw_info(ax_drs, bam1_label, u5o, u3o, orf_len, u5n, u3n, u5d, u3d,
                  len(d["drs_reads"]))
        ax_drs.set_xlim(xmin, xmax)
        ax_drs.set_xticks([])

        # SR reads
        draw_reads(ax_sr, d["sr_reads"], d["sr_rows"], u5o,
                   C_DRS2_FWD, C_DRS2_REV)
        draw_vlines(ax_sr, orf_len, u5n, u3n, u5d, u3d)
        draw_info(ax_sr, bam2_label, u5o, u3o, orf_len, u5n, u3n, u5d, u3d,
                  len(d["sr_reads"]))
        ax_sr.set_xlim(xmin, xmax)
        ax_sr.set_xlabel("Position relative to ORF start (bp)")

    # Shared legend (top right)
    handles = [
        mpatches.Patch(color=C_DRS1_FWD, label=f"{bam1_label} — fwd read"),
        mpatches.Patch(color=C_DRS1_REV, label=f"{bam1_label} — rev read"),
        mpatches.Patch(color=C_DRS2_FWD, label=f"{bam2_label} — fwd read"),
        mpatches.Patch(color=C_DRS2_REV, label=f"{bam2_label} — rev read"),
        mpatches.Patch(color=C_OUR,      label="Our annotation"),
        mpatches.Patch(color=C_NAG,      label="Reference annotation"),
        plt.Line2D([0], [0], color="black",   lw=1.2, ls="--",
                   label="Our ORF boundary"),
        plt.Line2D([0], [0], color=C_NAG,     lw=1.2, ls=":",
                   label="Reference boundary"),
        plt.Line2D([0], [0], color=C_DRS_BND, lw=1.2, ls=(0, (4, 2)),
                   label="DRS boundary (pre-merge)"),
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=FS_BODY,
               frameon=False, bbox_to_anchor=(0.99, 0.999),
               ncol=1, handlelength=1.6, handletextpad=0.5)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        path = outdir / f"{outname}.{ext}"
        fig.savefig(path, bbox_inches="tight",
                    dpi=300 if ext == "png" else None)
        print(f"Saved: {path}")
    plt.close(fig)


def main():
    args = parse_args()
    genes  = args.genes
    labels = args.panel_labels or [chr(ord("A") + i) for i in range(len(genes))]
    if len(labels) != len(genes):
        raise SystemExit("--panel_labels must match --genes in count")

    utr_our_df = load_utr_our(args.utr_our)
    utr_nag_df = load_utr_nag(args.utr_nag)
    utr_drs_df = load_utr_drs(args.utr_drs)

    gene_data = [
        prepare_gene(g, args.drs_bam, args.sr_bam,
                     utr_our_df, utr_nag_df, utr_drs_df,
                     args.name_suffix, args.max_reads, args.seed)
        for g in genes
    ]

    make_figure(genes, labels, gene_data, args.outdir, args.outname,
                bam1_label=args.bam1_label, bam2_label=args.bam2_label)


if __name__ == "__main__":
    main()
