#!/usr/bin/env python3
"""
plot_figure4.py
===============
Figure 4: Metagene coverage profiles around ORF boundaries.

Four panels (2 rows × 2 cols):
  Ai  SR  — short-read coverage anchored at ORF start (5' boundary)
  Aii SR  — short-read coverage anchored at ORF end   (3' boundary)
  Bi  DRS — DRS coverage anchored at ORF start
  Bii DRS — DRS coverage anchored at ORF end

Profiles are pre-computed TSV files with columns:
  position  mean_coverage  sem_coverage  n_transcripts

Expected filename pattern (in --profdir):
  {DT}_{ref}_{POINT}.tsv
  where DT  ∈ {SR_merged, DRS},
        ref ∈ {orf_only, orf_1000, our_ref, nagalakshmi},
        POINT ∈ {TSS, TES}

Edit REF_KEY_MAP below to match your file naming convention.

Outputs: figure4.svg  figure4.pdf  figure4.png

Usage:
    python plot_figure4.py \\
        --profdir /path/to/Metrics/profiles \\
        --outdir  /path/to/Figures
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# EDIT THESE TO MATCH YOUR DATA
# ══════════════════════════════════════════════════════════════════════════════

PALETTE = {
    "ORF only":    "#EE6677",
    "ORF ±1000":   "#CCBB44",
    "Our ref":     "#228833",
    "Nagalakshmi": "#AA3377",
}
REF_ORDER = ["ORF only", "ORF ±1000", "Our ref", "Nagalakshmi"]

# Map the reference fragment in the filename stem → display label
REF_KEY_MAP = {
    "orf_only":    "ORF only",
    "orf_1000":    "ORF ±1000",
    "our_ref":     "Our ref",
    "nagalakshmi": "Nagalakshmi",
}

DT_LABELS = {
    "SR":  "Short Read (Illumina)",
    "DRS": "DRS (Oxford Nanopore)",
}

# ══════════════════════════════════════════════════════════════════════════════
# TYPOGRAPHY
# ══════════════════════════════════════════════════════════════════════════════

FS_BODY   = 8
FS_LABEL  = 9
FS_PANEL  = 11
FS_TITLE  = 9
FS_LEGEND = 8

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FS_BODY,
    "axes.labelsize":    FS_LABEL,
    "axes.titlesize":    FS_TITLE,
    "xtick.labelsize":   FS_BODY,
    "ytick.labelsize":   FS_BODY,
    "axes.linewidth":    0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def parse_stem(stem):
    """Parse filename stem → (data_type, ref_label, point) or (None, None, None)."""
    for point in ("TSS", "TES"):
        if stem.endswith(f"_{point}"):
            name = stem[: -(len(point) + 1)]
            break
    else:
        return None, None, None

    if name.startswith("SR_merged_"):
        dt, rest = "SR",  name[len("SR_merged_"):]
    elif name.startswith("SR_"):
        dt, rest = "SR",  name[3:]
    elif name.startswith("DRS_"):
        dt, rest = "DRS", name[4:]
    else:
        return None, None, None

    label = REF_KEY_MAP.get(rest)
    return dt, label, point


def load_all(profdir):
    """Load all TSV profiles. Returns dict[(data_type, ref_label, point)] → DataFrame."""
    profiles = {}
    for tsv in sorted(Path(profdir).glob("*.tsv")):
        dt, label, point = parse_stem(tsv.stem)
        if label is None:
            continue
        df = pd.read_csv(tsv, sep="\t")
        df.columns = ["position", "mean_coverage", "sem_coverage", "n_transcripts"]
        if df.empty or df["mean_coverage"].isna().all():
            print(f"WARNING: empty profile {tsv.name}")
            continue
        profiles[(dt, label, point)] = df
        n = int(df["n_transcripts"].iloc[0])
        print(f"  {tsv.name}: n={n}, max={df['mean_coverage'].max():.1f}")
    return profiles

# ══════════════════════════════════════════════════════════════════════════════
# PANEL DRAWING
# ══════════════════════════════════════════════════════════════════════════════

def _panel_label(ax, label):
    ax.text(-0.14, 1.07, label, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", va="bottom", ha="left",
            clip_on=False)


def draw_metagene(ax, profiles, data_type, point,
                  panel_label, title, ylabel=None, show_legend=False):
    """Overlay metagene lines for all references on a single axes."""
    plotted = False
    for ref in REF_ORDER:
        key = (data_type, ref, point)
        if key not in profiles:
            continue
        df = profiles[key]
        ax.plot(df["position"], df["mean_coverage"],
                color=PALETTE[ref], lw=1.6, label=ref, zorder=3)
        plotted = True

    if not plotted:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        return

    ax.axvline(0, color="black", lw=1.0, ls="--", alpha=0.7, zorder=4)

    xlo, xhi = ax.get_xlim()
    if point == "TSS":
        ax.axvspan(0, xhi, alpha=0.06, color="steelblue", zorder=0)
        ax.text(40, ax.get_ylim()[1] * 0.97, "ORF \u2192",
                fontsize=FS_BODY - 1, color="steelblue", va="top", ha="left",
                clip_on=True)
    else:
        ax.axvspan(xlo, 0, alpha=0.06, color="steelblue", zorder=0)
        ax.text(-40, ax.get_ylim()[1] * 0.97, "\u2190 ORF",
                fontsize=FS_BODY - 1, color="steelblue", va="top", ha="right",
                clip_on=True)

    ax.set_xlabel("Distance from ORF boundary (bp)")
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=FS_TITLE, pad=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if show_legend:
        ax.legend(fontsize=FS_LEGEND, frameon=False,
                  loc="upper left", handlelength=1.4, handletextpad=0.4)

    _panel_label(ax, panel_label)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE ASSEMBLY
# ══════════════════════════════════════════════════════════════════════════════

def build_figure(profiles):
    fig = plt.figure(figsize=(13.0, 9.0), facecolor="white")
    gs  = fig.add_gridspec(2, 2,
                           hspace=0.50, wspace=0.32,
                           left=0.08, right=0.97,
                           top=0.95, bottom=0.10)

    ax_Ai  = fig.add_subplot(gs[0, 0])
    ax_Aii = fig.add_subplot(gs[0, 1])
    ax_Bi  = fig.add_subplot(gs[1, 0])
    ax_Bii = fig.add_subplot(gs[1, 1])

    # Shared y-limits within each data type
    def _ymax(dt):
        vals = [profiles[k]["mean_coverage"].max() for k in profiles if k[0] == dt]
        return max(vals) * 1.10 if vals else 1.0

    sr_ymax  = _ymax("SR")
    drs_ymax = _ymax("DRS")

    draw_metagene(ax_Ai,  profiles, "SR",  "TSS", "Ai",
                  "5\u2032 ORF boundary (ORF start)",
                  ylabel="Mean coverage (raw read depth)", show_legend=True)
    draw_metagene(ax_Aii, profiles, "SR",  "TES", "Aii",
                  "3\u2032 ORF boundary (ORF end)")
    draw_metagene(ax_Bi,  profiles, "DRS", "TSS", "Bi",
                  "5\u2032 ORF boundary (ORF start)",
                  ylabel="Mean coverage (raw read depth)", show_legend=True)
    draw_metagene(ax_Bii, profiles, "DRS", "TES", "Bii",
                  "3\u2032 ORF boundary (ORF end)")

    ax_Ai.set_ylim(0, sr_ymax);  ax_Aii.set_ylim(0, sr_ymax)
    ax_Bi.set_ylim(0, drs_ymax); ax_Bii.set_ylim(0, drs_ymax)

    # Row labels in left margin
    for ax, label in [(ax_Ai, DT_LABELS["SR"]), (ax_Bi, DT_LABELS["DRS"])]:
        ax.text(-0.16, 0.5, label, transform=ax.transAxes,
                fontsize=FS_BODY + 0.5, fontweight="bold",
                va="center", ha="right", rotation=90, clip_on=False)

    return fig

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profdir", required=True,
                   help="Directory containing *_TSS.tsv and *_TES.tsv profile files")
    p.add_argument("--outdir",  required=True,
                   help="Output directory for figure4.svg/pdf/png")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("Loading metagene profiles...")
    profiles = load_all(args.profdir)
    if not profiles:
        raise SystemExit("ERROR: no profiles loaded — check --profdir and filename pattern.")

    print("Building figure 4...")
    fig = build_figure(profiles)
    for ext, kw in [("svg", {}), ("pdf", {}), ("png", {"dpi": 300})]:
        out = outdir / f"figure4.{ext}"
        fig.savefig(out, bbox_inches="tight", **kw)
        print(f"Saved {out}")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
