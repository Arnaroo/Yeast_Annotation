#!/usr/bin/env python3
"""
plot_figure1a_schematic.py
=============================
Figure 1 panel A: a schematic of the segmentation pipeline (raw
coverage -> candidate breakpoints -> selected breakpoints -> final
annotation) for three synthetic genes. Nothing in it is measured and
the axes carry no numbers -- it illustrates the process, not a result;
the caption must say so.

Inputs: none, the panel is synthetic.
Output: figure1a.{pdf,png,svg}

Run:
  python3 plot_figure1a_schematic.py --outdir Figures_out
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

C_OUR = "#228833"; C_CDS = "#4477AA"; C_SENSE = "#333333"
C_CAND = "#BBBBBB"; C_SEL = "#EE7733"
FS_BODY, FS_LABEL, FS_PANEL = 7, 8, 10

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FS_BODY, "axes.labelsize": FS_LABEL, "xtick.labelsize": FS_BODY,
    "ytick.labelsize": FS_BODY, "legend.fontsize": FS_BODY, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none",
})

# Three genes in the schematic window. Each is a DRS coverage ramp: reads
# are captured at the poly(A) tail and extend toward the 5' end with
# variable length, so cumulative depth rises across the body and falls
# sharply at the 5' boundary.
GENES = [(0.06, 0.28, 0.62, "gene 3"), (0.36, 0.66, 1.00, "gene 2"), (0.74, 0.94, 0.48, "gene 1")]
BASELINE = 0.04
SELECTED = [0.06, 0.28, 0.36, 0.66, 0.74, 0.94]
SPURIOUS = [0.12, 0.17, 0.21, 0.25, 0.32, 0.41, 0.46, 0.50, 0.55, 0.60, 0.70, 0.79, 0.84, 0.89]
Y_TOP, Y_LIM, Y_NOTE = 1.10, 1.34, 1.22


def coverage(x):
    y = np.full_like(x, BASELINE)
    for start, end, peak, _ in GENES:
        inside = (x >= start) & (x <= end)
        frac = (x[inside] - start) / (end - start)
        y[inside] = BASELINE + (peak - BASELINE) * frac
    return y


def draw_coverage(ax):
    x = np.linspace(0.0, 1.0, 2000)
    ax.fill_between(x, 0.0, coverage(x), color=C_SENSE, lw=0.0)


def style(ax, title):
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, Y_LIM)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("Genomic position", labelpad=2)
    ax.set_title(title, fontsize=FS_BODY, pad=3)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fig = plt.figure(figsize=(7.0, 1.95))
    left, right, gap = 0.052, 0.995, 0.050
    w = (right - left - 3 * gap) / 4.0
    bottom, height = 0.175, 0.615
    xs = [left + i * (w + gap) for i in range(4)]
    axes = [fig.add_axes([x, bottom, w, height]) for x in xs]

    draw_coverage(axes[0])
    style(axes[0], "Input")

    draw_coverage(axes[1])
    n_cand = len(SELECTED) + len(SPURIOUS)
    axes[1].vlines(sorted(SELECTED + SPURIOUS), 0.0, Y_TOP, color=C_CAND, lw=0.6, zorder=3)
    style(axes[1], "1) Segmentation")
    axes[1].text(0.99, Y_NOTE, f"{n_cand} candidate breakpoints", ha="right", va="center",
                fontsize=FS_BODY - 1.0, color="#777777")

    draw_coverage(axes[2])
    axes[2].vlines(SPURIOUS, 0.0, Y_TOP, color=C_CAND, lw=0.4, alpha=0.35, zorder=2)
    axes[2].vlines(SELECTED, 0.0, Y_TOP, color=C_SEL, lw=1.1, zorder=4)
    style(axes[2], "2) Breakpoint selection")
    axes[2].text(0.99, Y_NOTE, f"{len(SELECTED)} selected breakpoints", ha="right", va="center",
                fontsize=FS_BODY - 1.0, color=C_SEL)

    draw_coverage(axes[3])
    for start, end, _peak, label in GENES:
        axes[3].fill_between([start, end], 0.0, Y_TOP, color=C_OUR, alpha=0.16, lw=0.0, zorder=1)
        axes[3].vlines([start, end], 0.0, Y_TOP, color=C_OUR, lw=0.9, zorder=4)
        axes[3].text((start + end) / 2.0, Y_TOP * 0.93, label, ha="center", va="center",
                    fontsize=FS_BODY - 1.0, color=C_CDS)
    style(axes[3], "3) Final annotation")
    axes[3].text(0.99, Y_NOTE, "expressed regions", ha="right", va="center",
                fontsize=FS_BODY - 1.0, color=C_OUR)

    fig.text(0.011, bottom + height / 2.0, "Number of reads", rotation=90, ha="center", va="center",
             fontsize=FS_LABEL)
    y_arrow = bottom + height / 2.0
    for i in range(3):
        x_from, x_to = xs[i] + w + 0.008, xs[i + 1] - 0.008
        fig.add_artist(mpatches.FancyArrowPatch((x_from, y_arrow), (x_to, y_arrow), transform=fig.transFigure,
                                                arrowstyle="-|>", mutation_scale=8, lw=0.8, color=C_SENSE,
                                                shrinkA=0, shrinkB=0))
    fig.text(0.004, 0.965, "A", fontsize=FS_PANEL, fontweight="bold", ha="left", va="top")

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(args.outdir, f"figure1a.{ext}"), dpi=400)
    plt.close(fig)
    print(f"wrote {args.outdir}/figure1a.{{pdf,png,svg}}")
    print("Panel A is a schematic. No axis carries a number and nothing in it is measured.")


if __name__ == "__main__":
    main()
