#!/usr/bin/env python3
"""
plot_figure2.py
===============
Figure 2: UTR annotation validation — scatter plots, back-to-back histogram,
and representative read pileup examples.

Layout
------
Row 0 (top):    Ai  Scatter: reference vs final 5\'  UTR length (Spearman r)
                Aii Scatter: reference vs final 3\'  UTR length (Spearman r)
                B   Back-to-back histogram: extension over reference (5\' left, 3\' right)
Row 1 (bottom): C, D, E  Read pileup examples (side-by-side)

All input paths are required CLI arguments — no hard-coded defaults.

Usage
-----
python plot_figure2.py \\
    --drs_utr    DRS_UTR_corrected.tsv \\
    --ref_utr    reference_UTR.csv \\
    --final_utr  final_utr.tsv \\
    --bam        mapped_reads.sorted.bam \\
    --genes      GENE_A GENE_B GENE_C \\
    --panel_labels C D E \\
    [--name_suffix _mRNA] \\
    [--outdir    ./figures] \\
    [--outname   figure2]

Input file formats
------------------
drs_utr   : TSV, columns: gene  max_five_prime_utr  max_three_prime_utr
ref_utr   : CSV, columns: Gene  five_prime_utr  three_prime_utr  (NA allowed)
final_utr : TSV, columns: Gene  final_five_prime_utr  final_three_prime_utr

Requirements: numpy, pandas, scipy, matplotlib>=3.4, pysam
"""

import argparse
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.transforms import blended_transform_factory
from scipy import stats
from pathlib import Path

# ── Colour palette ─────────────────────────────────────────────────────────────
C_FWD  = "#4477AA"   # forward-strand reads
C_REV  = "#EE6677"   # reverse-strand reads
C_OUR  = "#228833"   # our / DRS annotation (green)
C_REF  = "#AA3377"   # reference boundary retained (purple)
C_GREY = "#AAAAAA"   # within ±threshold (grey)

THR = 20             # nt threshold for DRS vs Naga categories

# ── Typography (applied globally via rcParams) ─────────────────────────────────
FS_BODY   = 8
FS_LABEL  = 9
FS_PANEL  = 11       # bold panel letter
FS_TITLE  = 9

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FS_BODY,
    "axes.labelsize":   FS_LABEL,
    "axes.titlesize":   FS_TITLE,
    "xtick.labelsize":  FS_BODY,
    "ytick.labelsize":  FS_BODY,
    "legend.fontsize":  FS_BODY,
    "axes.linewidth":   0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--drs_utr",     required=True,
                   help="Pre-merge DRS UTR annotation TSV")
    p.add_argument("--ref_utr",     required=True,
                   help="Reference UTR annotation CSV (NAs allowed)")
    p.add_argument("--final_utr",   required=True,
                   help="Final merged UTR annotation TSV")
    p.add_argument("--bam",         required=True,
                   help="Sorted, indexed BAM file")
    p.add_argument("--genes",       nargs="+", required=True,
                   help="Gene names for pileup panels")
    p.add_argument("--panel_labels", nargs="+", default=["C", "D", "E"])
    p.add_argument("--name_suffix", default="_mRNA")
    p.add_argument("--outdir",      default="./figures")
    p.add_argument("--outname",     default="figure2")
    p.add_argument("--max_reads",   type=int, default=300)
    p.add_argument("--seed",        type=int, default=0)
    return p.parse_args()


# ── Data loading ───────────────────────────────────────────────────────────────

def load_data(drs_file, ref_file, final_file):
    """Load and merge DRS, reference, and final UTR tables."""
    drs = pd.read_csv(drs_file, sep="\t")
    drs.columns = ["Gene", "drs5", "drs3"]

    nag = pd.read_csv(ref_file)
    nag.columns = ["Gene", "nag5", "nag3"]
    nag["nag5"] = pd.to_numeric(nag["nag5"], errors="coerce")
    nag["nag3"] = pd.to_numeric(nag["nag3"], errors="coerce")

    fin = pd.read_csv(final_file, sep="\t")
    fin.columns = ["Gene", "fin5", "fin3"]

    df = drs.merge(nag, on="Gene").merge(fin, on="Gene")
    print(f"Merged table: {len(df)} genes")
    return df


def make_scatter_df(df, end):
    """Genes with both Nagalakshmi and final UTR > 0 for scatter plots."""
    nag_col, fin_col = (f"nag{end}", f"fin{end}")
    sub = df[["Gene", nag_col, fin_col]].copy()
    sub.columns = ["Gene", "nag", "fin"]
    sub = sub.dropna(subset=["nag", "fin"])
    sub = sub[(sub["nag"] > 0) & (sub["fin"] > 0)]
    return sub.reset_index(drop=True)


def make_hist_df(df, end):
    """Genes with both final and Nagalakshmi UTR > 0; diff = final − Naga (≥ 0 by merge)."""
    cat_labels = [f"Within \u00b1{THR} nt", f"Extended (>{THR} nt)"]
    nag_col, fin_col = (f"nag{end}", f"fin{end}")
    sub = df[["Gene", nag_col, fin_col]].copy()
    sub.columns = ["Gene", "nag", "fin"]
    sub = sub.dropna(subset=["nag", "fin"])
    sub = sub[(sub["nag"] > 0) & (sub["fin"] > 0)].copy()
    sub["diff"] = sub["fin"] - sub["nag"]   # always ≥ 0 by merge construction
    sub["category"] = pd.cut(sub["diff"],
                             bins=[-np.inf, THR, np.inf],
                             labels=cat_labels)
    return sub.reset_index(drop=True), cat_labels

# ── Panel A: scatter plots (Nagalakshmi vs final) ──────────────────────────────

def _add_panel_label(ax, label):
    """Bold panel letter in upper-left corner, outside the axes."""
    ax.text(-0.14, 1.07, label, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", va="bottom", ha="left",
            clip_on=False)


def draw_scatter(ax, sdf, panel_label, title):
    """Scatter: x = Nagalakshmi UTR length, y = final UTR length (log scale)."""
    x = sdf["nag"].values
    y = sdf["fin"].values

    # Color points by whether DRS extended the boundary
    extended     = y > x
    not_extended = ~extended

    ax.scatter(x[not_extended], y[not_extended], c=C_REF, s=3,
               alpha=0.35, linewidths=0, rasterized=True,
               label="Reference boundary retained")
    ax.scatter(x[extended],     y[extended],     c=C_OUR, s=3,
               alpha=0.35, linewidths=0, rasterized=True,
               label="DRS extended boundary")

    # y = x identity line
    lo = 1.0
    hi = max(x.max(), y.max()) * 1.5
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.9, ls="--", zorder=5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    ax.set_xlabel("Reference UTR length (nt)")
    ax.set_ylabel("Final UTR length (nt)")
    ax.set_title(title, fontsize=FS_TITLE, pad=4)

    # Spearman r
    r, p = stats.spearmanr(x, y)
    pstr  = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    ax.text(0.05, 0.97, f"$r_s$ = {r:.2f}  ({pstr})",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=FS_BODY, style="italic",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=2))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _add_panel_label(ax, panel_label)


# ── Panel B: back-to-back histogram (5′ left | 3′ right) ──────────────────────

def draw_double_histogram(ax, hdf5, hdf3, cat_labels, panel_label, xlim_cap=600):
    """
    Butterfly / back-to-back histogram:
      Left  (x < 0): 5′ UTR extension, mirrored so bars grow leftward
      Right (x > 0): 3′ UTR extension
    y-axis: gene count (log scale); x-tick labels show absolute nt values.
    """
    from matplotlib.ticker import FuncFormatter

    cat_colours = {cat_labels[0]: C_GREY, cat_labels[1]: C_OUR}
    bins_pos = np.arange(0, xlim_cap + 10, 10)   # 0 … 600
    bins_neg = np.arange(-xlim_cap, 1, 10)        # -600 … 0

    # Draw grey (within ±THR) first, then green (extended) on top — both sides
    for cat in cat_labels:
        col = cat_colours[cat]
        vals5 = hdf5.loc[hdf5["category"] == cat, "diff"].clip(0, xlim_cap)
        vals3 = hdf3.loc[hdf3["category"] == cat, "diff"].clip(0, xlim_cap)
        ax.hist(-vals5, bins=bins_neg, color=col, alpha=0.85, edgecolor="none", zorder=2)
        ax.hist( vals3, bins=bins_pos, color=col, alpha=0.85, edgecolor="none", zorder=2)

    # Reference lines
    ax.axvline(0,    color="black",   lw=1.2, zorder=5)
    ax.axvline(-THR, color="#555555", lw=0.7, ls="--", zorder=4)
    ax.axvline( THR, color="#555555", lw=0.7, ls="--", zorder=4)

    ax.set_yscale("log")
    ax.set_xlim(-xlim_cap, xlim_cap)
    ax.set_xlabel("Extension over reference annotation (nt)")
    ax.set_ylabel("Number of genes")
    ax.set_title("5\u2032 UTR \u2190   |   \u2192 3\u2032 UTR", fontsize=FS_TITLE, pad=4)

    # Show absolute nt values on both sides of x-axis
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: str(int(abs(x)))))

    # Count annotations (blended: x in data, y in axes fraction)
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for hdf, x_sign, ha in [(hdf5, -1, "left"), (hdf3, 1, "right")]:
        n_ext    = int((hdf["category"] == cat_labels[1]).sum())
        n_within = int((hdf["category"] == cat_labels[0]).sum())
        ax.text(x_sign * xlim_cap * 0.97, 0.97,
                f"n = {n_ext}", transform=trans,
                va="top", ha=ha, fontsize=FS_BODY - 0.5,
                fontweight="bold", color=C_OUR)
        ax.text(x_sign * xlim_cap * 0.97, 0.52,
                f"n = {n_within}", transform=trans,
                va="top", ha=ha, fontsize=FS_BODY - 0.5,
                fontweight="bold", color="#777777")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _add_panel_label(ax, panel_label)

# ── Pileup helpers (adapted from plot_read_pileup_multigene.py) ────────────────

def _aligned_blocks(read):
    """Reference-coordinate aligned blocks, split at N (intron) CIGAR ops."""
    blocks, ref_pos = [], read.reference_start
    cur_start, cur_len = ref_pos, 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):         # M, =, X
            cur_len += length; ref_pos += length
        elif op == 2:               # D — deletion, stays in block
            cur_len += length; ref_pos += length
        elif op == 3:               # N — intron: break block
            blocks.append((cur_start, cur_start + cur_len))
            ref_pos += length; cur_start = ref_pos; cur_len = 0
        # I, S, H, P — not on reference, skip
    if cur_len > 0:
        blocks.append((cur_start, cur_start + cur_len))
    return blocks


def _get_ref_len(bam_path, contig):
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return dict(zip(bam.references, bam.lengths)).get(contig)


def _fetch_reads(bam_path, contig, ref_len, max_reads, seed):
    reads = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        try:
            for r in bam.fetch(contig, 0, ref_len):
                if r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                blocks = _aligned_blocks(r)
                if not blocks:
                    continue
                reads.append({"blocks": blocks,
                              "strand": "-" if r.is_reverse else "+",
                              "start":  blocks[0][0],
                              "end":    blocks[-1][1]})
        except ValueError:
            pass
    if len(reads) > max_reads:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(reads), size=max_reads, replace=False)
        reads = [reads[i] for i in sorted(idx)]
    return reads


def _assign_rows(reads, gap=15):
    reads_sorted = sorted(reads, key=lambda r: r["start"])
    row_ends = []
    for r in reads_sorted:
        placed = False
        for i, end in enumerate(row_ends):
            if r["start"] > end + gap:
                r["row"] = i; row_ends[i] = r["end"]; placed = True; break
        if not placed:
            r["row"] = len(row_ends); row_ends.append(r["end"])
    return reads_sorted, len(row_ends)


def prepare_gene_data(gene, bam_path, utr_final_df, name_suffix, max_reads, seed):
    """Fetch reads and compute geometry for one gene (no drawing)."""
    contig = gene + name_suffix
    if gene not in utr_final_df.index:
        return {"gene": gene, "error": f"{gene}: not in final UTR table"}
    utr5 = int(utr_final_df.loc[gene, "fin5"])
    utr3 = int(utr_final_df.loc[gene, "fin3"])
    ref_len = _get_ref_len(bam_path, contig)
    if ref_len is None:
        return {"gene": gene, "error": f"{gene}: contig '{contig}' not in BAM"}
    orf_len = ref_len - utr5 - utr3
    reads, n_rows = _assign_rows(_fetch_reads(bam_path, contig, ref_len, max_reads, seed))
    print(f"  {gene}: 5'UTR={utr5}  ORF={orf_len}  3'UTR={utr3}  "
          f"→ {len(reads)} reads, {n_rows} rows")
    return {"gene": gene, "error": None, "utr5": utr5, "utr3": utr3,
            "orf_len": orf_len, "reads": reads, "n_rows": n_rows}

def draw_pileup_panel(ax_annot, ax_reads, panel_label, data):
    """Draw annotation strip + read pileup for one gene."""
    gene = data["gene"]

    if data["error"]:
        ax_reads.text(0.5, 0.5, data["error"], ha="center", va="center",
                      transform=ax_reads.transAxes, color="grey", fontsize=FS_BODY)
        for ax in (ax_annot, ax_reads):
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_visible(False)
        ax_annot.text(-0.03, 0.5, panel_label, transform=ax_annot.transAxes,
                      fontsize=FS_PANEL, fontweight="bold", ha="right", va="center")
        return

    utr5, utr3, orf_len = data["utr5"], data["utr3"], data["orf_len"]
    reads, n_rows = data["reads"], data["n_rows"]
    READ_H = 0.8

    # — Annotation strip —
    if utr5 > 0:
        ax_annot.add_patch(mpatches.Rectangle((-utr5, 0), utr5, READ_H,
            facecolor=C_OUR, alpha=0.35, edgecolor=C_OUR, lw=0.8))
    ax_annot.add_patch(mpatches.Rectangle((0, 0), orf_len, READ_H,
        facecolor=C_OUR, alpha=0.85, edgecolor=C_OUR, lw=0.8))
    if utr3 > 0:
        ax_annot.add_patch(mpatches.Rectangle((orf_len, 0), utr3, READ_H,
            facecolor=C_OUR, alpha=0.35, edgecolor=C_OUR, lw=0.8))
    ax_annot.set_ylim(READ_H + 0.1, -0.1)
    ax_annot.set_yticks([])
    ax_annot.set_xticks([])
    for sp in ax_annot.spines.values(): sp.set_visible(False)

    # Panel label left of strip; gene name centred inside the green annotation bar
    ax_annot.text(-0.03, 0.5, panel_label, transform=ax_annot.transAxes,
                  fontsize=FS_PANEL, fontweight="bold", ha="right", va="center")
    x_center = (-utr5 + orf_len + utr3) / 2   # midpoint of transcript
    ax_annot.text(x_center, READ_H / 2, gene,
                  ha="center", va="center", fontsize=FS_BODY,
                  fontweight="bold", color="white", zorder=10)

    # — Read pileup —
    for r in reads:
        color = C_FWD if r["strand"] == "+" else C_REV
        y, x0, x1 = r["row"], r["blocks"][0][0] - utr5, r["blocks"][-1][1] - utr5
        ax_reads.plot([x0, x1], [y + READ_H / 2] * 2,
                      color=color, lw=0.5, alpha=0.45, zorder=1)
        for b0, b1 in r["blocks"]:
            ax_reads.add_patch(mpatches.Rectangle(
                (b0 - utr5, y), b1 - b0, READ_H,
                facecolor=color, edgecolor="none", alpha=0.70, zorder=2))
    ax_reads.set_ylim(-0.5, max(n_rows, 1) + 0.5)
    ax_reads.invert_yaxis()
    ax_reads.set_yticks([])

    # ORF region shading and boundary lines
    ax_reads.axvspan(0, orf_len, alpha=0.05, color="steelblue", zorder=0)
    ax_reads.axvline(0,       color="black", lw=0.8, ls="--", alpha=0.5)
    ax_reads.axvline(orf_len, color="black", lw=0.8, ls="--", alpha=0.5)

    # Info label
    ax_reads.text(0.01, 0.03,
                  f"5\u2032UTR = {utr5} \u00b7 ORF = {orf_len} \u00b7 3\u2032UTR = {utr3} nt"
                  f"   (n = {len(reads)} reads)",
                  transform=ax_reads.transAxes, fontsize=FS_BODY - 1,
                  fontweight="bold", va="bottom", ha="left",
                  bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))

    ax_reads.spines["top"].set_visible(False)
    ax_reads.spines["right"].set_visible(False)
    ax_reads.spines["left"].set_visible(False)

    ax_reads.set_xlabel("Position relative to ORF start (bp)", fontsize=FS_LABEL)

    xmin = -utr5 - 80
    xmax =  orf_len + utr3 + 80
    ax_annot.set_xlim(xmin, xmax)
    ax_reads.set_xlim(xmin, xmax)

# ── Figure assembly ────────────────────────────────────────────────────────────

def build_figure(df, gene_data_list, panel_labels_cde):
    """
    Figure layout (2 rows x 3 cols):
      Row 0: Ai scatter (5') | Aii scatter (3') | B back-to-back histogram
      Row 1: C pileup       | D pileup          | E pileup
    CDE legend sits in the hspace gap between the two rows.
    """
    FIG_W = 13.0; FIG_H = 11.5
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

    # 2-row × 3-col main gridspec
    gs = fig.add_gridspec(
        2, 3,
        height_ratios=[3.8, 6.5],
        hspace=0.42, wspace=0.38,
        left=0.09, right=0.97,
        top=0.97, bottom=0.07,
    )

    # ── Row 0: scatter panels and double histogram ────────────────────────────
    ax_Ai  = fig.add_subplot(gs[0, 0])
    ax_Aii = fig.add_subplot(gs[0, 1])
    ax_B   = fig.add_subplot(gs[0, 2])

    draw_scatter(ax_Ai,  make_scatter_df(df, "5"), "Ai",  "5\u2032 UTR")
    draw_scatter(ax_Aii, make_scatter_df(df, "3"), "Aii", "3\u2032 UTR")

    hdf5, cat_labels = make_hist_df(df, "5")
    hdf3, _          = make_hist_df(df, "3")
    draw_double_histogram(ax_B, hdf5, hdf3, cat_labels, "B")

    # Compact scatter/histogram legend embedded in the lower-right of ax_Aii
    # (lower right is empty because final ≥ Naga always → no points below y=x)
    scatter_handles = [
        mpatches.Patch(color=C_OUR, label="DRS extended"),
        mpatches.Patch(color=C_REF, label="Ref. boundary retained"),
        plt.Line2D([0], [0], color="black", lw=0.9, ls="--", label="y = x"),
        mpatches.Patch(color=C_GREY, label="Within \u00b120 nt"),
    ]
    ax_Aii.legend(handles=scatter_handles, loc="lower right",
                  frameon=False, fontsize=FS_BODY - 0.5, ncol=1,
                  handlelength=1.4, handletextpad=0.5)

    # ── Row 1: side-by-side pileup panels ────────────────────────────────────
    ANNOT_RATIO = 0.13   # annotation strip height as fraction of bottom row

    for col, (d, lbl) in enumerate(zip(gene_data_list, panel_labels_cde)):
        gs_sub = gs[1, col].subgridspec(2, 1,
                                          height_ratios=[ANNOT_RATIO, 1 - ANNOT_RATIO],
                                          hspace=0.0)
        ax_a = fig.add_subplot(gs_sub[0])
        ax_r = fig.add_subplot(gs_sub[1])
        draw_pileup_panel(ax_a, ax_r, lbl, d)

    # CDE legend in the hspace gap between the two rows (~figure y ≈ 0.63)
    cde_handles = [
        mpatches.Patch(color=C_FWD, label="Forward-strand read"),
        mpatches.Patch(color=C_REV, label="Reverse-strand read"),
        mpatches.Patch(color=C_OUR, label="Our annotation"),
    ]
    fig.legend(handles=cde_handles,
               loc="upper center",
               bbox_to_anchor=(0.5, 0.635),
               bbox_transform=fig.transFigure,
               ncol=3, frameon=False, fontsize=FS_BODY)

    return fig

# ── Main entry point ───────────────────────────────────────────────────────────

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes  = args.genes
    labels = args.panel_labels
    if len(labels) != len(genes):
        raise SystemExit("--panel_labels must match --genes in count")

    # ── Load annotation tables ────────────────────────────────────────────────
    print("Loading annotation tables...")
    df = load_data(args.drs_utr, args.ref_utr, args.final_utr)

    # Build a final-UTR indexed dataframe for pileup coordinate lookup.
    # Load directly from final_utr.tsv (not the merged df) so that genes absent
    # from Nagalakshmi are still available for pileup panels.
    fin_raw = pd.read_csv(args.final_utr, sep="\t")
    fin_raw.columns = ["Gene", "fin5", "fin3"]
    utr_final_indexed = fin_raw.set_index("Gene")

    # ── Prefetch reads for CDE genes (pass 1: geometry only) ─────────────────
    print("Fetching reads for pileup panels...")
    gene_data_list = [
        prepare_gene_data(g, args.bam, utr_final_indexed,
                          args.name_suffix, args.max_reads, args.seed)
        for g in genes
    ]

    # ── Build and save figure ─────────────────────────────────────────────────
    print("Building figure...")
    fig = build_figure(df, gene_data_list, labels)

    out_base = outdir / args.outname
    for ext, kw in [("svg",  {}),
                    ("pdf",  {}),
                    ("png",  {"dpi": 300})]:
        outpath = out_base.with_suffix(f".{ext}")
        fig.savefig(outpath, bbox_inches="tight", **kw)
        print(f"Saved: {outpath}")

    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
