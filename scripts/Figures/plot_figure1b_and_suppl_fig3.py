#!/usr/bin/env python3
"""
plot_figure1b_and_suppl_fig3.py
==================================
Figure 1 panel B and Supplementary Figure S3: the segmentation worked
example, redrawn on genomic coordinates with a coding-sequence track
and the Nagalakshmi annotation overlaid (Reviewer 1's specific request
for both figures).

Figure 1B follows one locus (default YBL091C) through four steps:
raw coverage (Bi), the full set of candidate breakpoints (Bii), the
breakpoints actually selected (Biii), and the resulting expressed
region against the final, post-merge boundary (Biv), with an
annotation track underneath showing this study's and Nagalakshmi's
UTRs.

Supplementary Figure S3 repeats the coverage + annotation pair for
eight further loci chosen to be difficult (a close neighbour, an
intron-containing gene), not representative.

Coverage source: sense = summed real per-position coverage across all
six construction conditions, in each gene's own transcript direction.
Antisense = the real sense coverage of any opposite-strand neighbour
gene whose own window overlaps this one, mapped into this window's
coordinates.

Input:  locus_coverage.tsv    panel, position, depth_sense, depth_antisense
        locus_features.tsv    panel, source(final/nagalakshmi), feature_gene,
                               feature_class(gene/CDS/intron/five_prime_UTR/
                               three_prime_UTR), strand, start, end
        locus_windows.tsv     panel, gene, chrom, strand, window_start, window_end
        locus_breakpoints.tsv panel, stage(candidate/selected), genomic_position
        locus_expressed.tsv   panel, expressed_start, expressed_end
        (all produced by locus_windows.py + locus_segmentation.R)
Output: figure1b.{pdf,png,svg}
        suppl_fig3.{pdf,png,svg}
        locus_utr_concordance.tsv

Run:
  python3 plot_figure1b_and_suppl_fig3.py --indir Data/LocusTracks --outdir Figures_out
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import blended_transform_factory

C_OUR = "#228833"; C_REF = "#AA3377"; C_CDS = "#4477AA"
C_SENSE = "#333333"; C_ANTI = "#EE6677"; C_CAND = "#BBBBBB"; C_SEL = "#EE7733"; C_GREY = "#AAAAAA"
FS_BODY, FS_LABEL, FS_PANEL = 7, 8, 10

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FS_BODY, "axes.labelsize": FS_LABEL, "xtick.labelsize": FS_BODY,
    "ytick.labelsize": FS_BODY, "legend.fontsize": FS_BODY, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none",
})

ROW_CDS, ROW_FIN, ROW_NAG = 2.0, 1.0, 0.0
BOX_H, UTR_H = 0.42, 0.30


def thousands(x, _pos=None):
    return f"{int(x):,}"


def draw_annotation(ax, feats, lo, hi, label_rows=True, window_gene=None, target_only_utrs=False):
    fin = feats[feats["source"] == "final"]
    nag = feats[feats["source"] == "nagalakshmi"]
    labels = []
    for gene, grp in fin.groupby("feature_gene"):
        cds = grp[grp["feature_class"] == "CDS"]
        if cds.empty:
            continue
        strand = cds["strand"].iloc[0]
        g_lo, g_hi = int(cds["start"].min()), int(cds["end"].max())
        ax.plot([g_lo, g_hi], [ROW_CDS, ROW_CDS], color=C_CDS, lw=0.7, zorder=2)
        for _, iv in grp[grp["feature_class"] == "intron"].iterrows():
            mid = (int(iv["start"]) + int(iv["end"])) / 2.0
            ax.plot(mid, ROW_CDS, marker=("<" if strand == "-" else ">"), ms=2.4, color=C_CDS, zorder=3)
        for _, ex in cds.iterrows():
            s, e = int(ex["start"]), int(ex["end"])
            ax.add_patch(mpatches.Rectangle((s, ROW_CDS - BOX_H / 2), e - s + 1, BOX_H,
                                            facecolor=C_CDS, edgecolor="none", zorder=4))
        tip = g_lo if strand == "-" else g_hi
        if lo <= tip <= hi:
            ax.plot(tip, ROW_CDS, marker=("<" if strand == "-" else ">"), ms=4.0, color=C_CDS, zorder=5)
        vis_lo, vis_hi = max(g_lo, lo), min(g_hi, hi)
        if vis_hi - vis_lo > (hi - lo) * 0.04:
            labels.append(((vis_lo + vis_hi - 2 * lo) / (2.0 * (hi - lo)), f"{gene} ({strand})", gene == window_gene))

    fs = FS_BODY - 1.2
    ax_pts = ax.get_position().width * ax.figure.get_figwidth() * 72.0
    xtrans = blended_transform_factory(ax.transAxes, ax.transData)
    labels.sort()
    for k, (xf, txt, is_win) in enumerate(labels):
        half = 0.5 * len(txt) * fs * 0.58 / ax_pts
        xf = min(max(xf, half), 1.0 - half)
        ax.text(xf, ROW_CDS + BOX_H / 2 + (0.16 if k % 2 == 0 else 0.62), txt, transform=xtrans,
                ha="center", va="bottom", fontsize=fs, color=C_CDS,
                fontweight=("bold" if is_win else "normal"), zorder=6)

    for table, row, colour in ((fin, ROW_FIN, C_OUR), (nag, ROW_NAG, C_REF)):
        guide = table[table["feature_class"] == "gene"]
        if target_only_utrs:
            guide = guide[guide["feature_gene"] == window_gene]
        for _, g in guide.iterrows():
            ax.plot([int(g["start"]), int(g["end"])], [row, row], color=C_GREY, lw=0.5, alpha=0.6, zorder=1)
        utrs = table[table["feature_class"].str.endswith("UTR")]
        if target_only_utrs:
            utrs = utrs[utrs["feature_gene"] == window_gene]
        for _, u in utrs.iterrows():
            s, e = int(u["start"]), int(u["end"])
            five = u["feature_class"] == "five_prime_UTR"
            ax.add_patch(mpatches.Rectangle((s, row - UTR_H / 2), max(e - s + 1, 1), UTR_H,
                                            facecolor=("white" if five else colour), edgecolor=colour,
                                            linewidth=0.55, zorder=4))
    if label_rows:
        suffix = f" ({window_gene})" if (target_only_utrs and window_gene) else ""
        for row, txt, colour in ((ROW_CDS, "CDS", C_CDS), (ROW_FIN, "This study" + suffix, C_OUR),
                                 (ROW_NAG, "Nagalakshmi" + suffix, C_REF)):
            ax.text(-0.008, row, txt, transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=FS_BODY - 0.5, color=colour)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.75, ROW_CDS + 1.30)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.xaxis.set_major_formatter(FuncFormatter(thousands))


def linear_strands(ax, sense_max, anti_max, pad=1.28):
    top, bot = max(int(sense_max), 1), max(int(anti_max), 1)
    m = max(top, bot) * pad
    ax.set_ylim(-m, m)
    ax.set_yticks([-bot, 0, top])
    ax.set_yticklabels([thousands(bot), "0", thousands(top)])


def style_cov(ax, lo, hi):
    ax.set_xlim(lo, hi)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0, labelbottom=False)


def panel_letter(ax, txt, dx=-0.095, dy=1.28):
    ax.text(dx, dy, txt, transform=ax.transAxes, ha="left", va="top", fontsize=FS_PANEL, fontweight="bold")


def load(indir):
    cov = pd.read_csv(os.path.join(indir, "locus_coverage.tsv"), sep="\t")
    feat = pd.read_csv(os.path.join(indir, "locus_features.tsv"), sep="\t")
    win = pd.read_csv(os.path.join(indir, "locus_windows.tsv"), sep="\t").set_index("panel")
    brk = pd.read_csv(os.path.join(indir, "locus_breakpoints.tsv"), sep="\t")
    exp = pd.read_csv(os.path.join(indir, "locus_expressed.tsv"), sep="\t")
    return cov, feat, win, brk, exp


def concordance_table(feat):
    u = feat[feat["feature_class"].str.endswith("UTR")]
    idx = ["panel", "window_gene", "feature_gene", "feature_class", "strand"]
    p = u.pivot_table(index=idx, columns="source", values=["start", "end"], aggfunc="first")
    p.columns = [f"{a}_{b}" for a, b in p.columns]
    p = p.reset_index()
    both = p.dropna(subset=["start_final", "start_nagalakshmi"])
    same = both[(both["start_final"] == both["start_nagalakshmi"]) & (both["end_final"] == both["end_nagalakshmi"])]
    p["status"] = np.where(p["start_nagalakshmi"].isna(), "final_only",
                           np.where(p["start_final"].isna(), "nagalakshmi_only",
                                   np.where((p["start_final"] == p["start_nagalakshmi"]) &
                                            (p["end_final"] == p["end_nagalakshmi"]), "identical", "different")))
    print(f"  UTRs called by both sources : {len(both)}")
    print(f"  identical coordinates       : {len(same)}")
    print(f"  called only by this study   : {int((p['status'] == 'final_only').sum())}")
    print(f"  called only by Nagalakshmi  : {int((p['status'] == 'nagalakshmi_only').sum())}")
    return p.sort_values(["panel", "feature_gene", "feature_class"])


def coverage_arrays(cov, panel):
    s = cov[cov["panel"] == panel].sort_values("position")
    return s["position"].to_numpy(), s["depth_sense"].to_numpy(), s["depth_antisense"].to_numpy()


def draw_fig1b(cov, feat, win, brk, exp, outdir, panel_name):
    w = win.loc[panel_name]
    lo, hi = int(w["window_start"]), int(w["window_end"])
    pos, sense, anti = coverage_arrays(cov, panel_name)
    fe = feat[feat["panel"] == panel_name]
    cand = brk[(brk["panel"] == panel_name) & (brk["stage"] == "candidate")]
    sel = brk[(brk["panel"] == panel_name) & (brk["stage"] == "selected")]
    ex = exp[exp["panel"] == panel_name]

    fig = plt.figure(figsize=(7.0, 6.2))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.25, 1, 1, 1, 1.05], hspace=0.36,
                          left=0.14, right=0.985, top=0.945, bottom=0.115)
    ax_i, ax_ii, ax_iii, ax_iv = [fig.add_subplot(gs[k]) for k in range(4)]
    ax_ann = fig.add_subplot(gs[4])

    ax_i.fill_between(pos, 0, sense, color=C_SENSE, lw=0, zorder=3)
    ax_i.fill_between(pos, 0, -anti, color=C_ANTI, lw=0, alpha=0.85, zorder=2)
    ax_i.axhline(0, color="black", lw=0.6, zorder=4)
    linear_strands(ax_i, sense.max(), anti.max())
    ax_i.set_ylabel("Reads\n(both strands)")
    style_cov(ax_i, lo, hi)
    panel_letter(ax_i, "Bi")
    ax_i.legend(handles=[mpatches.Patch(color=C_SENSE, label=f"sense ({w['strand']} strand)"),
                         mpatches.Patch(color=C_ANTI, alpha=0.85, label="antisense")],
                loc="lower right", frameon=False, ncol=2, handlelength=1.1, borderpad=0.1,
                columnspacing=1.0, bbox_to_anchor=(1.0, 1.0))

    for ax, letter in ((ax_ii, "Bii"), (ax_iii, "Biii"), (ax_iv, "Biv")):
        ax.fill_between(pos, 0, sense, color=C_SENSE, lw=0, zorder=3)
        ax.set_ylim(0, sense.max() * 1.28)
        style_cov(ax, lo, hi)
        panel_letter(ax, letter)
    box_top, box_bot = ax_ii.get_position().y1, ax_iv.get_position().y0
    fig.text(0.028, (box_top + box_bot) / 2.0, "Reads (sense strand)", rotation=90,
             ha="center", va="center", fontsize=FS_LABEL)

    edge_positions = set()
    for _, r in ex.iterrows():
        edge_positions.add(int(r["expressed_start"]))
        edge_positions.add(int(r["expressed_end"]))
    cand_pos = cand["genomic_position"].tolist()
    edge_cands = [x for x in cand_pos if x in edge_positions]

    for x in cand_pos:
        is_edge = x in edge_positions
        ax_ii.axvline(x, color=(C_OUR if is_edge else C_CAND), lw=(1.3 if is_edge else 0.5),
                     zorder=(4 if is_edge else 1))
    ax_ii.text(0.995, 0.93, f"{len(cand)} candidate breakpoints", transform=ax_ii.transAxes,
              ha="right", va="top", fontsize=FS_BODY - 0.5, color="#777777")
    if edge_cands:
        ax_ii.text(0.995, 0.80, "boundary reached here (see Biv)", transform=ax_ii.transAxes,
                  ha="right", va="top", fontsize=FS_BODY - 0.5, color=C_OUR)

    for x in sel["genomic_position"]:
        is_edge = x in edge_positions
        ax_iii.axvline(x, color=(C_OUR if is_edge else C_SEL), lw=(1.3 if is_edge else 0.9), zorder=4)
    ax_iii.text(0.995, 0.93, f"{len(sel)} selected breakpoints", transform=ax_iii.transAxes,
               ha="right", va="top", fontsize=FS_BODY - 0.5, color=C_SEL)
    not_in_selected = sorted(set(edge_cands) - set(sel["genomic_position"]))
    if not_in_selected:
        ax_iii.text(0.995, 0.80, "some Biv edges are not selected here (reached via Bii instead)",
                   transform=ax_iii.transAxes, ha="right", va="top", fontsize=FS_BODY - 0.5, color=C_GREY)

    for _, r in ex.iterrows():
        ax_iv.axvspan(r["expressed_start"], r["expressed_end"], color=C_OUR, alpha=0.18, lw=0, zorder=1)
        ax_iv.axvline(r["expressed_start"], color=C_OUR, lw=0.9, zorder=4)
        ax_iv.axvline(r["expressed_end"], color=C_OUR, lw=0.9, zorder=4)
    ax_iv.text(0.995, 0.93, "DRS call (pre-merge)", transform=ax_iv.transAxes, ha="right", va="top",
              fontsize=FS_BODY - 0.5, color=C_OUR)

    for x in edge_cands:
        con = mpatches.ConnectionPatch(xyA=(x, ax_ii.get_ylim()[0]), coordsA=ax_ii.transData,
                                       xyB=(x, ax_iv.get_ylim()[1]), coordsB=ax_iv.transData,
                                       color=C_OUR, lw=0.6, ls=":", alpha=0.6, zorder=0)
        fig.add_artist(con)

    wgene_final = fe[(fe["source"] == "final") & (fe["feature_gene"] == w["gene"]) &
                     (fe["feature_class"].isin(["CDS", "five_prime_UTR", "three_prime_UTR"]))]
    if not wgene_final.empty:
        f_lo, f_hi = int(wgene_final["start"].min()), int(wgene_final["end"].max())
        for x in (f_lo, f_hi):
            ax_iv.axvline(x, color=C_OUR, lw=1.1, ls="--", zorder=5)
        ax_iv.text(0.995, 0.78, "final boundary (post-merge)", transform=ax_iv.transAxes,
                  ha="right", va="top", fontsize=FS_BODY - 0.5, color=C_OUR, fontstyle="italic")

    draw_annotation(ax_ann, fe, lo, hi, window_gene=w["gene"], target_only_utrs=True)
    ax_ann.set_xlabel(f"Chromosome {w['chrom']} position (nt)")
    fig.legend(handles=[mpatches.Patch(facecolor=C_CDS, edgecolor="none", label="CDS"),
                        mpatches.Patch(facecolor="white", edgecolor="black", label="5' UTR"),
                        mpatches.Patch(facecolor="black", edgecolor="black", label="3' UTR")],
              loc="lower center", ncol=3, frameon=False, fontsize=FS_BODY, handlelength=1.1,
              columnspacing=1.6, bbox_to_anchor=(0.5, 0.004))

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(outdir, f"figure1b.{ext}"), dpi=400)
    plt.close(fig)
    print(f"{panel_name} {w['gene']} {w['chrom']}:{lo}-{hi}  candidate {len(cand)}  "
          f"selected {len(sel)}  expressed {len(ex)}")


def draw_suppl_fig3(cov, feat, win, exp, outdir, prefix="S1"):
    panels = [p for p in win.index if p.startswith(prefix)]
    ncol, nrow = 2, (len(panels) + 1) // 2
    fig = plt.figure(figsize=(7.2, 9.4))
    outer = fig.add_gridspec(nrow, ncol, hspace=0.46, wspace=0.20, left=0.115, right=0.99,
                             top=0.962, bottom=0.085)
    fig.legend(handles=[
        mpatches.Patch(color=C_SENSE, label="sense coverage"),
        mpatches.Patch(color=C_ANTI, alpha=0.85, label="antisense coverage"),
        mpatches.Patch(color=C_OUR, alpha=0.25, label="DRS call (pre-merge)"),
        plt.Line2D([0], [0], color=C_OUR, lw=1.1, ls="--", label="final boundary (post-merge)"),
        mpatches.Patch(facecolor="white", edgecolor="black", label="5' UTR"),
        mpatches.Patch(facecolor="black", edgecolor="black", label="3' UTR"),
    ], loc="lower center", ncol=6, frameon=False, handlelength=1.4, columnspacing=1.3,
        bbox_to_anchor=(0.5, 0.002))

    for k, panel in enumerate(panels):
        r, c = divmod(k, ncol)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[1.6, 1.0], hspace=0.06)
        ax_c = fig.add_subplot(inner[0])
        ax_a = fig.add_subplot(inner[1])
        w = win.loc[panel]
        lo, hi = int(w["window_start"]), int(w["window_end"])
        pos, sense, anti = coverage_arrays(cov, panel)
        ex = exp[exp["panel"] == panel]
        pfeat = feat[feat["panel"] == panel]

        for _, rr in ex.iterrows():
            ax_c.axvspan(rr["expressed_start"], rr["expressed_end"], color=C_OUR, alpha=0.15, lw=0, zorder=1)
        ax_c.fill_between(pos, 0, sense, color=C_SENSE, lw=0, zorder=3)
        ax_c.fill_between(pos, 0, -anti, color=C_ANTI, lw=0, alpha=0.85, zorder=2)
        ax_c.axhline(0, color="black", lw=0.5, zorder=4)

        wgene_final = pfeat[(pfeat["source"] == "final") & (pfeat["feature_gene"] == w["gene"]) &
                            (pfeat["feature_class"].isin(["CDS", "five_prime_UTR", "three_prime_UTR"]))]
        if not wgene_final.empty:
            f_lo, f_hi = int(wgene_final["start"].min()), int(wgene_final["end"].max())
            for x in (f_lo, f_hi):
                ax_c.axvline(x, color=C_OUR, lw=1.1, ls="--", zorder=5)

        linear_strands(ax_c, sense.max(), anti.max())
        style_cov(ax_c, lo, hi)
        panel_letter(ax_c, panel[-1], dx=-0.10, dy=1.16)
        ax_c.set_title(f"{w['gene']} ({w['strand']})", fontsize=FS_BODY, pad=2.0)

        draw_annotation(ax_a, pfeat, lo, hi, label_rows=(c == 0), window_gene=w["gene"], target_only_utrs=True)
        ax_a.set_xlabel(f"Chr {w['chrom']} position (nt)", labelpad=1.5)
        ax_a.tick_params(axis="x", labelsize=FS_BODY - 1.5)
        for lbl in ax_a.get_xticklabels():
            lbl.set_rotation(20)
            lbl.set_ha("right")
        print(f"{panel:5s} {w['gene']:11s} {w['chrom']}:{lo}-{hi}  "
              f"final boundary drawn: {not wgene_final.empty}  expressed {len(ex)}")

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(outdir, f"suppl_fig3.{ext}"), dpi=400)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--fig1b-panel", default="Fig1B", help="panel id in locus_windows.tsv for Figure 1B")
    ap.add_argument("--suppl-prefix", default="S1", help="panel-id prefix for the eight-loci figure")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cov, feat, win, brk, exp = load(args.indir)
    print("UTR concordance across the drawn windows")
    conc = concordance_table(feat)
    conc.to_csv(os.path.join(args.outdir, "locus_utr_concordance.tsv"), sep="\t", index=False)

    draw_fig1b(cov, feat, win, brk, exp, args.outdir, args.fig1b_panel)
    print()
    draw_suppl_fig3(cov, feat, win, exp, args.outdir, args.suppl_prefix)
    print("\nwrote figure1b.{pdf,png,svg}")
    print("wrote suppl_fig3.{pdf,png,svg}")


if __name__ == "__main__":
    main()
