#!/usr/bin/env python3
"""
sensitivity_curves.py
========================
Supplementary Figure S1. One line per threshold, showing how much the
post-merge, six-condition maximum boundary call (the released quantity)
changes as that one threshold is moved across its full grid range with
the other three held at the released value. Complements the one-step
bar chart in sensitivity_per_parameter_onestep.tsv (which only compares
the released combination to its immediate neighbours) by showing the
shape of the response across all five grid points per parameter, not
just the first step.

Input:  sensitivity_whole_grid.tsv (sensitivity_analyse.py)
Output: sensitivity_curves.{pdf,png,svg}

Run:  python3 sensitivity_curves.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "..", "..", "Data", "Sensitivity")

REL = dict(min_reads=20, merge_window=40, effect_frac=0.10, expr_frac=0.05)
GRID = dict(
    min_reads=[5, 10, 20, 40, 80],
    merge_window=[10, 20, 40, 80, 160],
    effect_frac=[0.025, 0.05, 0.10, 0.20, 0.35],
    expr_frac=[0.01, 0.025, 0.05, 0.10, 0.20],
)
LABELS = dict(min_reads="coverage floor (reads)", merge_window="clustering window (nt)",
              effect_frac="effect size (frac of max)", expr_frac="expression floor (frac)")
COLORS = dict(min_reads="#4c72b0", merge_window="#55a868",
              effect_frac="#c44e52", expr_frac="#8172b2")


def one_parameter_slice(conc, param):
    others = [p for p in REL if p != param]
    mask = True
    for p in others:
        mask = mask & (conc[p] == REL[p])
    sub = conc.loc[mask].sort_values(param)
    assert len(sub) == 5, f"{param}: expected 5 rows, got {len(sub)}"
    return sub


def main():
    conc = pd.read_csv(os.path.join(DATA_DIR, "sensitivity_whole_grid.tsv"), sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.8), sharey=True)
    for end, ax, col in [("five", axes[0], "ident_five"), ("three", axes[1], "ident_three")]:
        for param in GRID:
            sub = one_parameter_slice(conc, param)
            x = list(range(-2, 3))  # grid step relative to the released value
            y = (1 - sub[col].values) * 100  # % of genes whose boundary changed
            ax.plot(x, y, marker="o", color=COLORS[param], label=LABELS[param], linewidth=2)
        ax.axvline(0, color="grey", linestyle=":", linewidth=1)
        ax.set_xticks(range(-2, 3))
        ax.set_xlabel("grid step from released value")
        ax.set_title(f"{end} end")
    axes[0].set_ylabel("% of genes with a changed boundary\n(vs released six-condition max)")
    axes[1].legend(loc="upper left", fontsize=9, framealpha=0.9)
    fig.suptitle("Boundary sensitivity across the full grid, one threshold varied at a time")
    fig.tight_layout()

    # Annotate actual grid values below each axis, since step units differ by parameter.
    lines = []
    for param in GRID:
        vals = ", ".join(str(v) for v in GRID[param])
        lines.append(f"{LABELS[param]}: {vals} (released {REL[param]})")
    fig.text(0.02, -0.12, "\n".join(lines), fontsize=13, ha="left", va="top")

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(DATA_DIR, f"sensitivity_curves.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote sensitivity_curves.{pdf,png,svg}")


if __name__ == "__main__":
    main()
