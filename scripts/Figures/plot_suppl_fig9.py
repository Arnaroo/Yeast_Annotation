#!/usr/bin/env python3
"""
plot_suppl_fig9.py
====================
Supplementary Figure S9. Metagene coverage split by whether this study
moved the boundary: the same our_ref genes and profiles as Figure 4,
split on whether the final annotation moved that boundary outward by
more than 20 nt relative to Nagalakshmi. Two rows per technology (mean
depth, then per-gene scaled depth), two columns (ORF start, ORF end).

Three drawing rules, all consequences of what the data can support:
  1. Curves are truncated where fewer than 20% of a group's genes still
     have reference sequence.
  2. Mean depth ignores positions where the sequence does not exist,
     rather than counting them as zero.
  3. The scaled statistic plotted is the median across genes, not the
     mean, since a mean of per-gene ratios is dominated by the smallest
     denominators.

Axis labels say ORF start and ORF end, never transcript start and end.

Input:  metagene_scaled_profiles.tsv (compute_metagene_scaled.py)
Output: suppl_fig9.{pdf,png,svg}
        metagene_scaled_contrasts.tsv

Run:
  python3 plot_suppl_fig9.py --profiles Data/Metrics/profiles_scaled/metagene_scaled_profiles.tsv --outdir Figures_out
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_EXT = "#228833"     # boundary extended, green
C_UNC = "#AA3377"     # boundary unchanged, purple
C_ORF = "#4477AA"

FS_BODY, FS_LABEL, FS_PANEL = 7, 8, 10
SUPPORT = 0.20  # minimum fraction of a group still having sequence

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FS_BODY,
    "axes.labelsize": FS_LABEL,
    "xtick.labelsize": FS_BODY,
    "ytick.labelsize": FS_BODY,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "svg.fonttype": "none",
})

ANCHOR_LABEL = {"TSS": "ORF start", "TES": "ORF end"}
DT_LABEL = {"SR": "Short read (Illumina)", "DRS": "DRS (Nanopore)"}
STAT_LABEL = {
    "depth": "Mean depth (reads)",
    "scaled_median": "Median scaled depth\n(gene ORF mean = 1)",
}


def supported(sub):
    n = sub["n_genes"].iloc[0]
    return (sub["n_at_position"] / max(n, 1)) >= SUPPORT


def draw_panel(ax, prof, dt, anchor, stat, show_legend=False):
    xlo, xhi = np.inf, -np.inf
    for group, colour, label in (("unchanged", C_UNC, "boundary unchanged"),
                                 ("extended", C_EXT, "boundary extended")):
        sub = prof[(prof.data_type == dt) & (prof.anchor == anchor) &
                   (prof.statistic == stat) & (prof.group == group)]
        sub = sub.sort_values("position")
        m = supported(sub).to_numpy()
        if not m.any():
            continue
        x = sub["position"].to_numpy()[m]
        y = sub["value"].to_numpy()[m]
        s = sub["spread"].to_numpy()[m]
        n = int(sub["n_genes"].iloc[0])
        ax.fill_between(x, y - s, y + s, color=colour, alpha=0.16, lw=0, zorder=2)
        ax.plot(x, y, color=colour, lw=1.2, zorder=3, label=f"{label} (n = {n:,})")
        xlo, xhi = min(xlo, x.min()), max(xhi, x.max())

    pad = 0.02 * (xhi - xlo)
    ax.set_xlim(xlo - pad, xhi + pad)
    ax.axvline(0, color="black", lw=0.7, ls="--", alpha=0.8, zorder=4)
    lo, hi = ax.get_xlim()
    if anchor == "TSS":
        ax.axvspan(0, hi, color=C_ORF, alpha=0.05, lw=0, zorder=0)
    else:
        ax.axvspan(lo, 0, color=C_ORF, alpha=0.05, lw=0, zorder=0)
    ax.set_xlabel(f"Distance from {ANCHOR_LABEL[anchor]} (nt)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(y=0.12)
    if show_legend:
        ax.legend(loc="upper left", frameon=False, fontsize=FS_BODY - 0.5,
                  handlelength=1.3, borderpad=0.1, labelspacing=0.25)


def panel_letter(ax, txt):
    ax.text(-0.20, 1.06, txt, transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", ha="left", va="bottom")


def contrasts(prof):
    rows = []
    for dt in ("SR", "DRS"):
        for anchor in ("TSS", "TES"):
            for stat in ("depth", "scaled_median"):
                sel = prof[(prof.data_type == dt) & (prof.anchor == anchor) &
                           (prof.statistic == stat)]
                piv = sel.pivot_table(index="position", columns="group", values="value")
                cnt = sel.pivot_table(index="position", columns="group", values="n_at_position")
                ng = sel.groupby("group")["n_genes"].first()
                for pos in (-300, -100, 100, 300):
                    if pos not in piv.index:
                        continue
                    e, u = piv.loc[pos, "extended"], piv.loc[pos, "unchanged"]
                    fe = cnt.loc[pos, "extended"] / ng["extended"]
                    fu = cnt.loc[pos, "unchanged"] / ng["unchanged"]
                    rows.append({
                        "data_type": dt, "anchor": ANCHOR_LABEL[anchor], "statistic": stat,
                        "position": pos, "extended": round(e, 4), "unchanged": round(u, 4),
                        "ratio": round(e / u, 3) if u else np.nan,
                        "n_extended_at_position": int(cnt.loc[pos, "extended"]),
                        "n_unchanged_at_position": int(cnt.loc[pos, "unchanged"]),
                        "supported": bool(fe >= SUPPORT and fu >= SUPPORT),
                    })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    prof = pd.read_csv(args.profiles, sep="\t")

    fig = plt.figure(figsize=(7.2, 8.4))
    gs = fig.add_gridspec(4, 2, hspace=0.62, wspace=0.30,
                          left=0.135, right=0.985, top=0.925, bottom=0.075)

    letters = ["Ai", "Aii", "Bi", "Bii", "Ci", "Cii", "Di", "Dii"]
    k = 0
    for r, (dt, stat) in enumerate([("SR", "depth"), ("SR", "scaled_median"),
                                    ("DRS", "depth"), ("DRS", "scaled_median")]):
        for c, anchor in enumerate(("TSS", "TES")):
            ax = fig.add_subplot(gs[r, c])
            draw_panel(ax, prof, dt, anchor, stat, show_legend=(r == 0 and c == 0))
            panel_letter(ax, letters[k])
            k += 1
            if c == 0:
                ax.set_ylabel(STAT_LABEL[stat])
            if r in (0, 2) and c == 0:
                ax.text(-0.30, 1.19, DT_LABEL[dt], transform=ax.transAxes,
                        fontsize=FS_LABEL, fontweight="bold", ha="left", va="bottom")

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(args.outdir, f"suppl_fig9.{ext}"), dpi=400)
    plt.close(fig)
    print(f"wrote {args.outdir}/suppl_fig9.{{pdf,png,svg}}")

    ct = contrasts(prof)
    ct_path = os.path.join(args.outdir, "metagene_scaled_contrasts.tsv")
    ct.to_csv(ct_path, sep="\t", index=False)
    print(f"wrote {ct_path} ({len(ct)} rows)")
    print(ct[ct.supported].to_string(index=False))


if __name__ == "__main__":
    main()
