#!/usr/bin/env python3
"""
plot_figure3.py
===============
Figure 3: Alignment quality comparison across reference transcriptomes.

Self-contained — no dependency on companion scripts.

Layout (2 rows × 4 panels, spacer between blocks):
  Row 0: Ai  SR mapping rate  | Aii DRS mapping rate  || Ci  SR antisense   | Cii DRS antisense
  Row 1: Bi  SR multi-mapping | Bii DRS multi-mapping || Di  DRS clip total  | Dii DRS clip frac

Required inputs (--metrics directory structure):
  flagstat/
    SR_fastq_read_counts.txt          — TSV: reference, total_input_reads, mapping_rate_pct
    {stem}.flagstat.txt               — samtools flagstat output
    {stem}.multimapping.txt           — TSV: total_primary, multi_mapped, multi_rate
    {stem}.bwa_multimapping.txt       — TSV: total_primary, mapq0_reads, mapq0_rate, method
  antisense/
    {stem}.antisense.txt              — TSV: primary_mapped, antisense
  aln_per_read/
    {stem}.aln_per_read.txt           — TSV: n_alignments, read_count
  softclip/
    {stem}.softclip.tsv               — TSV per DRS transcript: mean_clip_total, mean_clip_frac, …

Edit LABEL_MAP below to match your BAM/file stem naming convention.
STAR_TRANSCRIPTOME lists stems for which antisense is suppressed by STAR.
BWA_GENOME_STEMS lists bwa-mem2 genome stems whose multi-mapping proxy is MAPQ=0.

Outputs: figure3.svg  figure3.pdf  figure3.png

Usage:
    python plot_figure3.py \\
        --metrics /path/to/Metrics \\
        --outdir  /path/to/Figures
"""

import argparse, re
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
    "Genome":      "#4477AA",
    "ORF only":    "#EE6677",
    "ORF ±1000":   "#CCBB44",
    "Our ref":     "#228833",
    "Nagalakshmi": "#AA3377",
}
REF_ORDER = ["Genome", "ORF only", "ORF ±1000", "Our ref", "Nagalakshmi"]

# Maps file stem → (data_type, reference_label)
LABEL_MAP = {
    "DRS_genome":            ("DRS", "Genome"),
    "DRS_orf_only":          ("DRS", "ORF only"),
    "DRS_orf_1000":          ("DRS", "ORF ±1000"),
    "DRS_our_ref":           ("DRS", "Our ref"),
    "DRS_nagalakshmi":       ("DRS", "Nagalakshmi"),
    "SR_merged_genome":      ("SR",  "Genome"),
    "SR_merged_orf_only":    ("SR",  "ORF only"),
    "SR_merged_orf_1000":    ("SR",  "ORF ±1000"),
    "SR_merged_our_ref":     ("SR",  "Our ref"),
    "SR_merged_nagalakshmi": ("SR",  "Nagalakshmi"),
}

# STAR transcriptome BAMs suppress antisense — mark those stems here
STAR_TRANSCRIPTOME = {
    "SR_merged_orf_only", "SR_merged_orf_1000",
    "SR_merged_our_ref",  "SR_merged_nagalakshmi",
}

# bwa-mem2 genome BAMs: never write secondary records; use MAPQ=0 as proxy
BWA_GENOME_STEMS = {"SR_merged_genome", "DRS_genome"}

DT_LABELS = {"SR": "Short Read (Illumina)", "DRS": "DRS (Oxford Nanopore)"}

# ══════════════════════════════════════════════════════════════════════════════
# TYPOGRAPHY
# ══════════════════════════════════════════════════════════════════════════════

FS_BODY  = 8
FS_LABEL = 9
FS_PANEL = 11
FS_TITLE = 9
FS_ANNOT = 7.5

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

XLABELS = ["Genome", "ORF only", "ORF\n\u00b11000", "Our ref", "Nagalakshmi"]

# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def _parse_flagstat(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"(\d+) \+ \d+ (.+)", line.strip())
            if not m:
                continue
            val, key = int(m.group(1)), m.group(2).strip()
            if "in total" in key:      d["total_in_bam"]    = val
            elif key == "primary":     d["primary"]         = val
            elif key == "secondary":   d["secondary"]       = val
            elif key.startswith("primary mapped"):
                                       d["primary_mapped"]  = val
    if d.get("primary") and d.get("primary_mapped") is not None:
        d["mapping_rate_bam_pct"] = round(d["primary_mapped"] / d["primary"] * 100, 2)
    return d


def _load_multimapping(path):
    if not Path(path).exists():
        return {}
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if len(lines) < 2:
        return {}
    d = dict(zip(lines[0].split("\t"), lines[-1].split("\t")))
    if "multi_rate" in d:
        try:    return {"multi_rate_pct": round(float(d["multi_rate"]) * 100, 2)}
        except: return {}
    if "secondary_per_primary" in d:
        try:    return {"secondary_per_primary_pct": round(float(d["secondary_per_primary"]) * 100, 2)}
        except: return {}
    if "mapq0_rate" in d:
        try:    return {"mapq0_rate_pct": round(float(d["mapq0_rate"]) * 100, 2)}
        except: return {}
    return {}


def _load_antisense(path, stem):
    if stem in STAR_TRANSCRIPTOME:
        return {"antisense_rate_pct": float("nan")}
    if not Path(path).exists():
        return {}
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if len(lines) < 2:
        return {}
    d = dict(zip(lines[0].split("\t"), lines[-1].split("\t")))
    try:
        total = int(d["primary_mapped"])
        anti  = int(d["antisense"])
        return {"antisense_rate_pct": round(anti / total * 100, 2) if total else float("nan")}
    except (KeyError, ValueError):
        return {}


def _load_aln_per_read(path):
    if not Path(path).exists():
        return pd.DataFrame(columns=["n_alignments", "read_count"])
    try:
        df = pd.read_csv(path, sep="\t")
        df.columns = ["n_alignments", "read_count"]
        return df.dropna().astype(int).sort_values("n_alignments").reset_index(drop=True)
    except Exception:
        return pd.DataFrame(columns=["n_alignments", "read_count"])


def load_all_metrics(metrics_dir):
    mdir = Path(metrics_dir)

    # SR mapping rates from FASTQ counts file (more accurate denominator)
    sr_fastq = {}
    fastq_path = mdir / "flagstat" / "SR_fastq_read_counts.txt"
    if fastq_path.exists():
        tmp = pd.read_csv(fastq_path, sep="\t")
        for _, row in tmp.iterrows():
            sr_fastq[row["reference"]] = {
                "total_input":      int(row["total_input_reads"]),
                "mapping_rate_pct": float(row["mapping_rate_pct"]),
            }

    rows = []
    for stem, (dt, ref) in LABEL_MAP.items():
        row = {"data_type": dt, "reference": ref, "stem": stem}

        fp = mdir / "flagstat" / f"{stem}.flagstat.txt"
        if fp.exists():
            row.update(_parse_flagstat(fp))

        # Mapping rate
        if dt == "SR" and stem in sr_fastq:
            row["mapping_rate_pct"] = sr_fastq[stem]["mapping_rate_pct"]
            row["total_input"]      = sr_fastq[stem]["total_input"]
        elif dt == "DRS":
            row["mapping_rate_pct"] = row.get("mapping_rate_bam_pct")

        # Multi-mapping (NH-tag for STAR; MAPQ=0 proxy loaded separately)
        row.update(_load_multimapping(mdir / "flagstat" / f"{stem}.multimapping.txt"))
        row.update(_load_multimapping(mdir / "flagstat" / f"{stem}.bwa_multimapping.txt"))

        # For bwa-mem2 genome BAMs promote MAPQ=0 proxy to multi_rate_pct
        if stem in BWA_GENOME_STEMS and row.get("mapq0_rate_pct") is not None:
            row["multi_rate_pct"] = row["mapq0_rate_pct"]

        row.update(_load_antisense(mdir / "antisense" / f"{stem}.antisense.txt", stem))
        rows.append(row)

    return pd.DataFrame(rows)


def load_softclip_all(metrics_dir):
    mdir = Path(metrics_dir) / "softclip"
    dfs  = []
    for stem, (dt, ref) in LABEL_MAP.items():
        if dt != "DRS":
            continue
        fp = mdir / f"{stem}.softclip.tsv"
        if fp.exists():
            df = pd.read_csv(fp, sep="\t")
            df["reference"] = ref
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def enrich_drs_multimapping(df, metrics_dir):
    """Fill DRS multi_rate_pct from alignments-per-read distributions."""
    df   = df.copy()
    mdir = Path(metrics_dir) / "aln_per_read"
    for stem, (dt, ref) in LABEL_MAP.items():
        if dt != "DRS":
            continue
        aln = _load_aln_per_read(mdir / f"{stem}.aln_per_read.txt")
        if aln.empty:
            continue
        total  = aln["read_count"].sum()
        unique = aln.loc[aln["n_alignments"] == 1, "read_count"].sum()
        pct    = round((1 - unique / total) * 100, 2) if total else float("nan")
        mask   = (df["data_type"] == "DRS") & (df["reference"] == ref)
        df.loc[mask, "multi_rate_pct"] = pct
    return df

# ══════════════════════════════════════════════════════════════════════════════
# PANEL DRAWING
# ══════════════════════════════════════════════════════════════════════════════

def _panel_label(ax, label):
    ax.text(-0.16, 1.08, label, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", va="bottom", ha="left",
            clip_on=False)


def _group_header(ax, text):
    ax.text(0.0, 1.30, text, transform=ax.transAxes,
            fontsize=FS_BODY + 0.5, fontweight="bold",
            va="bottom", ha="left", clip_on=False)


def _spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_bar(ax, df, data_type, metric, panel_label, ylabel=None,
             title=None, ylim=None, na_label="N/A", genome_star=False):
    sub  = (df[df["data_type"] == data_type]
              .set_index("reference").reindex(REF_ORDER))
    vals = (sub[metric] if metric in sub.columns
            else pd.Series([float("nan")] * len(REF_ORDER), index=REF_ORDER))
    ymax = ylim[1] if ylim else (vals.dropna().max() * 1.25 if vals.dropna().any() else 1)

    bars = ax.bar(range(len(REF_ORDER)), vals.fillna(0),
                  color=[PALETTE[r] for r in REF_ORDER],
                  width=0.62, edgecolor="white", linewidth=0.8)

    for i, (bar, val) in enumerate(zip(bars, vals)):
        if pd.isna(val):
            ax.text(bar.get_x() + bar.get_width() / 2, ymax * 0.04,
                    na_label, ha="center", va="bottom",
                    fontsize=FS_ANNOT, color="grey", style="italic")
        else:
            txt = f"{val:.1f}" + ("*" if genome_star and REF_ORDER[i] == "Genome" else "")
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + ymax * 0.018,
                    txt, ha="center", va="bottom", fontsize=FS_ANNOT)

    if metric == "mapping_rate_pct":
        ax.axhline(100, color="grey", lw=0.8, ls="--", alpha=0.5)

    ax.set_xticks(range(len(REF_ORDER)))
    ax.set_xticklabels(XLABELS, rotation=35, ha="right", fontsize=FS_BODY)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(ylabel or "")
    if not ylabel:
        ax.set_yticklabels([])
    if title:
        ax.set_title(title, fontsize=FS_TITLE, pad=4)
    _spines(ax)
    _panel_label(ax, panel_label)


def draw_violin(ax, sc_df, metric, panel_label, ylabel=None, title=None):
    refs   = [r for r in REF_ORDER if r in sc_df["reference"].unique()]
    data   = [sc_df.loc[sc_df["reference"] == r, metric].dropna().values for r in refs]
    parts  = ax.violinplot(data, positions=range(len(refs)),
                           showmedians=True, showextrema=False)
    for pc, col in zip(parts["bodies"], [PALETTE[r] for r in refs]):
        pc.set_facecolor(col); pc.set_alpha(0.75); pc.set_edgecolor("white")
    parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.8)

    ax.figure.canvas.draw()
    ylo, yhi = ax.get_ylim(); yrange = yhi - ylo
    fmt = "{:.0f} bp" if metric == "mean_clip_total" else "{:.3f}"
    for i, d in enumerate(data):
        if len(d):
            ax.text(i, np.median(d) + yrange * 0.04, fmt.format(np.median(d)),
                    ha="center", va="bottom", fontsize=FS_ANNOT, fontweight="bold")

    ax.set_xticks(range(len(refs)))
    ax.set_xticklabels(refs, rotation=35, ha="right", fontsize=FS_BODY)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(ylabel or "")
    if not ylabel:
        ax.set_yticklabels([])
    if title:
        ax.set_title(title, fontsize=FS_TITLE, pad=4)
    _spines(ax)
    _panel_label(ax, panel_label)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE ASSEMBLY
# ══════════════════════════════════════════════════════════════════════════════

def build_figure(df, sc_df):
    def _ymax(metric, factor=1.25):
        v = df[metric].dropna(); return (v.max() * factor) if len(v) else 1.0

    map_ylim   = (0, 115)
    multi_ylim = (0, _ymax("multi_rate_pct",    1.25))
    anti_ylim  = (0, _ymax("antisense_rate_pct", 1.30))

    fig = plt.figure(figsize=(14.0, 9.5), facecolor="white")
    gs  = fig.add_gridspec(2, 5,
                           width_ratios=[1, 1, 0.10, 1, 1],
                           hspace=0.78, wspace=0.35,
                           left=0.07, right=0.97,
                           top=0.93, bottom=0.14)

    ax_Ai  = fig.add_subplot(gs[0, 0]); ax_Aii = fig.add_subplot(gs[0, 1])
    ax_Ci  = fig.add_subplot(gs[0, 3]); ax_Cii = fig.add_subplot(gs[0, 4])
    ax_Bi  = fig.add_subplot(gs[1, 0]); ax_Bii = fig.add_subplot(gs[1, 1])
    ax_Di  = fig.add_subplot(gs[1, 3]); ax_Dii = fig.add_subplot(gs[1, 4])
    for row in range(2):
        fig.add_subplot(gs[row, 2]).set_visible(False)

    _group_header(ax_Ai, "Mapping rate")
    _group_header(ax_Bi, "Multi-mapping rate\u2020")
    _group_header(ax_Ci, "Antisense alignment rate\u2021")
    _group_header(ax_Di, "Soft-clip distribution (DRS only)")

    draw_bar(ax_Ai,  df, "SR",  "mapping_rate_pct",   "Ai",
             ylabel="Reads mapped (%)", title="Short Read (Illumina)", ylim=map_ylim)
    draw_bar(ax_Aii, df, "DRS", "mapping_rate_pct",   "Aii",
             title="DRS (Oxford Nanopore)", ylim=map_ylim)
    draw_bar(ax_Bi,  df, "SR",  "multi_rate_pct",     "Bi",
             ylabel="Multi-mapping reads (%)", title="Short Read (Illumina)",
             ylim=multi_ylim, genome_star=True)
    draw_bar(ax_Bii, df, "DRS", "multi_rate_pct",     "Bii",
             title="DRS (Oxford Nanopore)", ylim=multi_ylim)
    draw_bar(ax_Ci,  df, "SR",  "antisense_rate_pct", "Ci",
             ylabel="Antisense reads (%)", title="Short Read (Illumina)",
             ylim=anti_ylim, na_label="N/A\n(STAR)")
    draw_bar(ax_Cii, df, "DRS", "antisense_rate_pct", "Cii",
             title="DRS (Oxford Nanopore)", ylim=anti_ylim)
    draw_violin(ax_Di,  sc_df, "mean_clip_total", "Di",
                ylabel="Mean total soft-clip per read (bp)", title="DRS")
    draw_violin(ax_Dii, sc_df, "mean_clip_frac",  "Dii",
                ylabel="Mean clip fraction", title="DRS")

    fig.text(0.02, 0.005,
             "\u2020 SR Genome: MAPQ=0 primary reads (* bar) as multi-mapping proxy "
             "(bwa-mem2 never writes 0x100 secondary records); "
             "SR transcriptomes: NH \u2265 2 (STAR); DRS: secondary alignments (minimap2).\n"
             "\u2021 STAR transcriptome BAMs: N/A \u2014 STAR suppresses antisense by default.",
             fontsize=FS_BODY - 1.5, color="#555555", va="bottom", linespacing=1.6)
    return fig

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metrics", required=True,
                   help="Root metrics directory (containing flagstat/, antisense/, etc.)")
    p.add_argument("--outdir",  required=True,
                   help="Output directory for figure3.svg/pdf/png")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("Loading metrics...")
    df    = load_all_metrics(args.metrics)
    sc_df = load_softclip_all(args.metrics)
    df    = enrich_drs_multimapping(df, args.metrics)

    print("Building figure 3...")
    fig = build_figure(df, sc_df)
    for ext, kw in [("svg", {}), ("pdf", {}), ("png", {"dpi": 300})]:
        out = outdir / f"figure3.{ext}"
        fig.savefig(out, bbox_inches="tight", **kw)
        print(f"Saved {out}")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
