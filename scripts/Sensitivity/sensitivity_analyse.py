#!/usr/bin/env python3
"""
sensitivity_analyse.py
=======================
Turns the raw per-condition sweep (sensitivity_sweep.R) into the two
tables the manuscript reports: the whole-grid summary (how much the
released, six-condition-max boundary call moves across all 625
threshold combinations) and the per-parameter one-step breakdown (which
of the four thresholds actually drives the movement).

The released DRS-derived boundary (DRS_UTR_corrected.tsv) is the
per-gene maximum across the six conditions (NS_PS, NS_RDT, S10_PS,
S10_RDT, NS, S10), so every comparison here is made on that
six-condition maximum at each grid point, not on any single condition's
own call -- a condition's own call can be insensitive to a threshold
move while the max across conditions still changes, if a different
condition becomes the new maximum.

Two views are produced:
  1. Whole-grid summary (one row per combination, 625 rows): for every
     combination, the fraction of genes whose call is retained,
     identical, or within 10 nt of the released combination's call.
     This is the primary robustness claim -- how far can a joint,
     worst-case move across all four thresholds shift the calls.
  2. Per-parameter one-step breakdown (one row per end x parameter):
     holding three thresholds at the released value and moving the
     fourth by a single grid step, how much of the gene population
     moves and by how much. The median move is 0 nt for every
     parameter and would misleadingly read as "no parameter matters"
     on its own; frac_moved and frac_moved_gt10nt show the minority
     that actually does.

Input:  sensitivity_sweep_summary.tsv, sensitivity_sweep_per_gene.tsv.gz
        (both from sensitivity_sweep.R)
Output: sensitivity_whole_grid.tsv
        sensitivity_per_parameter_onestep.tsv

Run:  python3 sensitivity_analyse.py
"""

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "..", "..", "Data", "Sensitivity")

REL = dict(min_reads=20, merge_window=40, effect_frac=0.10, expr_frac=0.05)
PARAM_COLS = ["min_reads", "merge_window", "effect_frac", "expr_frac"]


def log(msg=""):
    print(msg)


def main():
    ssum = pd.read_csv(os.path.join(DATA_DIR, "sensitivity_sweep_summary.tsv"), sep="\t")
    sgene = pd.read_csv(os.path.join(DATA_DIR, "sensitivity_sweep_per_gene.tsv.gz"), sep="\t")

    combos = ssum.drop_duplicates("combo")[["combo"] + PARAM_COLS]
    rel_mask = combos.apply(lambda r: all(np.isclose(r[c], REL[c]) for c in PARAM_COLS), axis=1)
    assert rel_mask.sum() == 1, "expected exactly one combination matching the released thresholds"
    rel_combo = int(combos.loc[rel_mask, "combo"].iloc[0])
    log("released settings are combination %d" % rel_combo)

    # The real, released quantity: six-condition max per gene, per combination.
    maxpg = sgene.groupby(["combo", "gene"])[["five", "three"]].max().reset_index()
    rel_calls = maxpg[maxpg.combo == rel_combo].set_index("gene")
    log("genes with a six-condition max call at the released settings: %d" % len(rel_calls))
    log("")

    # -- 1. Whole-grid summary -------------------------------------------
    rows = []
    for k in sorted(combos["combo"].unique()):
        d = maxpg[maxpg.combo == k].set_index("gene")
        common = rel_calls.index.intersection(d.index)
        if len(common) == 0:
            continue
        d5 = (d.loc[common, "five"] - rel_calls.loc[common, "five"]).abs()
        d3 = (d.loc[common, "three"] - rel_calls.loc[common, "three"]).abs()
        rows.append(dict(
            combo=k,
            retained=len(common) / len(rel_calls),
            identical=float(((d5 == 0) & (d3 == 0)).mean()),
            ident_five=float((d5 == 0).mean()),
            ident_three=float((d3 == 0).mean()),
            within10=float(((d5 <= 10) & (d3 <= 10)).mean()),
            n_genes=len(d),
            mean_five=d["five"].mean(),
            mean_three=d["three"].mean(),
        ))
    whole = combos.merge(pd.DataFrame(rows), on="combo").sort_values("combo")
    whole.to_csv(os.path.join(DATA_DIR, "sensitivity_whole_grid.tsv"), sep="\t", index=False)
    log("wrote sensitivity_whole_grid.tsv, %d rows" % len(whole))

    other = whole[whole["combo"] != rel_combo]
    log("  across the %d non-released combinations:" % len(other))
    log("    retained   %.3f to %.3f" % (other["retained"].min(), other["retained"].max()))
    log("    identical  %.3f to %.3f" % (other["identical"].min(), other["identical"].max()))
    log("    within10   %.3f to %.3f" % (other["within10"].min(), other["within10"].max()))
    log("")

    # -- 2. Per-parameter one-step breakdown ------------------------------
    def is_one_step(row):
        return sum(1 for c in PARAM_COLS if not np.isclose(row[c], REL[c])) == 1

    neighbours = combos[combos.apply(is_one_step, axis=1)]
    log("one-step neighbours of the released combination: %d" % len(neighbours))

    onestep_rows = []
    for _, nrow in neighbours.iterrows():
        k = int(nrow["combo"])
        moved_param = next(c for c in PARAM_COLS if not np.isclose(nrow[c], REL[c]))
        d = maxpg[maxpg.combo == k].set_index("gene")
        common = rel_calls.index.intersection(d.index)
        for end in ("five", "three"):
            diff = (d.loc[common, end] - rel_calls.loc[common, end]).abs()
            onestep_rows.append(pd.DataFrame({
                "combo": k, "moved_param": moved_param, "moved_to": nrow[moved_param],
                "end": end, "gene": common, "abs_diff_nt": diff.values,
            }))
    onestep = pd.concat(onestep_rows, ignore_index=True)

    summary_rows = []
    for end in ("five", "three"):
        for param in PARAM_COLS:
            d = onestep[(onestep.end == end) & (onestep.moved_param == param)]["abs_diff_nt"]
            if len(d) == 0:
                continue
            summary_rows.append(dict(
                end=end, parameter=param, n_genes=len(d),
                frac_moved=float((d > 0).mean()),
                frac_moved_gt10nt=float((d > 10).mean()),
                median_move_nt=float(d.median()),
                max_move_nt=float(d.max()),
            ))
    per_param = pd.DataFrame(summary_rows)
    per_param.to_csv(os.path.join(DATA_DIR, "sensitivity_per_parameter_onestep.tsv"),
                      sep="\t", index=False)
    log("wrote sensitivity_per_parameter_onestep.tsv, %d rows" % len(per_param))
    log("")
    log(per_param.to_string(index=False))
    log("")
    log("READING: a one-step nudge on any single one of the four reviewed thresholds "
        "leaves the great majority of genes' released (six-condition max) boundary "
        "exactly unchanged -- but not all: a real minority do move, some "
        "substantially. The median of 0 nt is a genuine majority effect, not "
        "evidence that the thresholds never matter; frac_moved and "
        "frac_moved_gt10nt are the numbers that show the minority.")


if __name__ == "__main__":
    main()
