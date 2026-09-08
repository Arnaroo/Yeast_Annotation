# A refined *Saccharomyces cerevisiae* reference transcriptome from Direct RNA Sequencing

This repository contains the data, code, and pre-built reference files associated with:

> **Rossini O, Cleynen A, Shirokikh N.** A refined *Saccharomyces cerevisiae* reference transcriptome from Direct RNA Sequencing, with a reusable pipeline for UTR annotation updates. *FEMS Yeast Research*, in revision (manuscript FEMSYR-26-06-0101).

The annotation, the transcript FASTA and GFF3, and the merged UTR table are also archived on Zenodo under DOI [10.5281/zenodo.20828039](https://doi.org/10.5281/zenodo.20828039).

The resource provides per-gene 5′ and 3′ UTR lengths for *S. cerevisiae* derived from Oxford Nanopore Direct RNA Sequencing (DRS), merged with the widely used Nagalakshmi *et al.* (2008) short-read RNA-seq annotation. The merge rule guarantees that **no Nagalakshmi boundary is ever shortened**: at each boundary and for each gene, the longer of the two estimates is retained. Genes detected only in DRS receive DRS-only values; genes absent from DRS retain their Nagalakshmi boundaries.

---

## Table of contents

- [Repository structure](#repository-structure)
- [Data files](#data-files)
- [Pipeline overview](#pipeline-overview)
  - [Step 1 — Aggregate DRS profiles (R language)](#step-1--aggregate-drs-profiles-r)
  - [Step 2 — Change-point segmentation (R language)](#step-2--change-point-segmentation-r)
  - [Step 3 — UTR table and merge (R language)](#step-3--utr-table-and-merge-r)
  - [Step 4 — Build transcript FASTA reference (Python)](#step-4--build-transcript-fasta-reference-python)
  - [Preprocessing — alignment and metrics (HPC shell scripts)](#preprocessing--alignment-and-metrics-hpc-shell-scripts)
  - [Figures](#figures)
  - [Comparisons — validation against independent annotations and references](#comparisons--validation-against-independent-annotations-and-references)
- [Sensitivity — segmentation threshold sweep](#sensitivity--segmentation-threshold-sweep)
- [Dependencies](#dependencies)
- [Quick start](#quick-start)
- [Citation](#citation)
- [License](#license)

---

## Repository structure

```
.
├── Data/
│   ├── Annotation/          # Final annotation tables and FASTA references
│   ├── Metrics/             # Pre-computed alignment quality metrics
│   │   ├── aln_per_read/
│   │   ├── antisense/
│   │   ├── flagstat/
│   │   ├── profiles/        # Metagene coverage profiles (TSV)
│   │   ├── profiles_scaled/ # compute_metagene_scaled.py / plot_suppl_fig9.py outputs
│   │   ├── softclip/
│   │   └── archive_pre_real_rerun/   # metrics as they stood before the real-rerun
│   ├── UNAGI/               # unagi_compare.py outputs
│   ├── Overlap/             # neighbour_overlap.py outputs
│   ├── Abundance/           # abundance_concordance.py / read_migration_analyse.py outputs
│   ├── TIFseq/              # tifseq_rates.py outputs
│   ├── Shortfall/           # shortfall.py / read_geometry.py outputs
│   ├── LocusTracks/         # locus_*.py/.R and Figure 1 / Supplementary Fig. S3 outputs
│   └── Sensitivity/         # sensitivity_sweep.R / sensitivity_analyse.py outputs
├── FinalReference/          # Ready-to-use FASTA references
├── scripts/
│   ├── ReferenceConstruction/   # Steps 1–4: build the annotation
│   │   ├── 01_aggregate_data.R
│   │   ├── 02_transcripts_segmentation.R
│   │   ├── 03_utr_table.R
│   │   └── 04_build_utr_reference.py
│   ├── Preprocessing/           # HPC scripts: alignment and metrics
│   │   ├── 01_build_STAR_indices.sh
│   │   ├── 02_align_DRS.sh
│   │   ├── 03_align_ShortRead.sh
│   │   ├── 04_metagene_coverage.sh
│   │   ├── 05_metrics.sh
│   │   ├── compute_metagene.py
│   │   └── softclip_metrics.py
│   ├── Figures/                 # Figure plotting scripts
│   │   ├── plot_figure2.py
│   │   ├── plot_figure3.py
│   │   ├── plot_figure4.py
│   │   ├── plot_figure1a_schematic.py
│   │   ├── locus_windows.py
│   │   ├── locus_coverage_real.py
│   │   ├── locus_segmentation.R
│   │   ├── plot_figure1b_and_suppl_fig3.py
│   │   ├── compose_figure1.py
│   │   ├── plot_suppl_fig2.py         # superseded, see the Figures table
│   │   ├── plot_suppl_fig3.R
│   │   ├── plot_suppl_fig4.py
│   │   ├── compute_metagene_scaled.py
│   │   └── plot_suppl_fig9.py
│   ├── Comparisons/             # Validation against independent annotations and references
│   │   ├── derive_pelechano_utr.py
│   │   ├── plot_final_vs_tifseq_comparison.R
│   │   ├── tifseq_rates.py
│   │   ├── plot_suppl_fig6.R
│   │   ├── shortfall.py
│   │   ├── read_geometry.py
│   │   ├── plot_suppl_fig7.R
│   │   ├── unagi_compare.py
│   │   ├── neighbour_overlap.py
│   │   ├── abundance_concordance.py
│   │   ├── read_migration.sh
│   │   └── read_migration_analyse.py
│   └── Sensitivity/              # Segmentation threshold sweep
│       ├── segmentation_corrected.R
│       ├── sensitivity_sweep.R
│       ├── sensitivity_analyse.py
│       └── sensitivity_curves.py
└── README.md
```

Note: `Data/Bam/` (example BAM files) was removed from the repository and from its history — the files were large enough to make cloning impractical and are better hosted alongside the other large data on Zenodo. `.gitignore` keeps it from being re-added by accident.

---

## Data files

### Core annotation tables (`Data/Annotation/`)

| File | Description |
|---|---|
| `DRS_UTR_corrected.tsv` | Pre-merge DRS annotation. Per-gene 5′ and 3′ UTR lengths (nt) derived from DRS, taking the maximum across all conditions and replicates. Columns: `gene`, `max_five_prime_utr`, `max_three_prime_utr`. |
| `Nagalakshmi_UTR.csv` | Nagalakshmi *et al.* (2008) short-read RNA-seq UTR lengths. Columns: `Gene`, `five_prime_utr`, `three_prime_utr`. |
| `final_utr.tsv` | **Final merged annotation** (recommended for most uses). Per-boundary maximum of `DRS_UTR_corrected.tsv` and `Nagalakshmi_UTR.csv`. Columns: `Gene`, `final_five_prime_utr`, `final_three_prime_utr`. |
| `gene_models_backbone_with_UTR_slots.gff3` | Curated gene-model backbone GFF3 (gene + CDS + intron features with UTR coordinate slots). Used by `04_build_utr_reference.py` to reconstruct the transcript FASTA. |
| `S288C_reference_sequence_R64-4-1_20230823.fsa` | SGD S288C genomic FASTA (release R64-4-1, 2023-08-23). |
| `Final_UTRs.gff3` | GFF3 with final UTR coordinate annotations overlaid on the gene-model backbone. |
| `Final_genome_chr_id.fa` / `.fai` | Genomic FASTA with renamed chromosome identifiers (I–XVI, mt). |
| `Final_w.fa` | **Final transcript FASTA reference** built from `final_utr.tsv` via `04_build_utr_reference.py`. One spliced mRNA per gene, including UTR sequence. |
| `Nagalakshmi_genome_chr_id.fa` / `.fai` | Genomic FASTA (chromosome IDs renamed) for the Nagalakshmi reference build. |
| `Nagalakshmi_UTRs.gff3` | GFF3 with Nagalakshmi UTR coordinates for the alternative reference build. |
| `Nagalakshmi_w.fa` | Transcript FASTA built from `Nagalakshmi_UTR.csv` (for comparison). |
| `TARGET_final_reference_w.fa` | Concatenated multi-reference FASTA used for competitive alignment benchmarking. |
| `Pelechano_UTR.csv` | Per-gene 5′/3′ UTR lengths derived from the Pelechano, Wei & Steinmetz (2013) TIF-seq major-isoform table, via `scripts/Comparisons/derive_pelechano_utr.py`. Used as an independent, orthogonal benchmark that took no part in building the annotation. Columns: `Gene`, `five_prime_utr`, `three_prime_utr`. |

### Metrics (`Data/Metrics/`)

Pre-computed alignment quality metrics for all five references × two data types, ready for direct use with the figure scripts:

- `flagstat/` — `samtools flagstat` output and multimapping rates
- `antisense/` — antisense alignment rates
- `aln_per_read/` — alignments-per-read distributions
- `softclip/` — soft-clip statistics (DRS only)
- `profiles/` — metagene coverage profiles around the ORF start and end (files keep their internal `_TSS`/`_TES` naming; see the note under `compute_metagene.py`)

---

## Pipeline overview

### Step 1 — Aggregate DRS profiles (R language)

**Script:** `scripts/ReferenceConstruction/01_aggregate_data.R`

Sums per-replicate DRS read-count profiles (one count per gene × nucleotide position) across all replicates within each condition, producing a single `{condition}_combined.tsv` per condition. Edit `DATA_DIR` and `CONDITIONS` before running.

**Input:** per-replicate TSV files `{condition}_{replicate}.txt` (columns: `gene`, `position`, `transcript_np`)
**Output:** `{condition}_combined.tsv`
**R packages:** `dplyr`, `readr`

---

### Step 2 — Change-point segmentation (R language)

**Script:** `scripts/ReferenceConstruction/02_transcripts_segmentation.R`

For each gene passing a minimum read-count threshold (default 20), runs change-point detection on the cumulative coverage profile using **Segmentor3IsBack** (negative-binomial model). Nearby breakpoints are merged and low-amplitude transitions discarded. The remaining breakpoints delimit "expressed" (high-coverage) regions, from which UTR boundaries are inferred.

Run once per condition. Edit `DATA_DIR`, `RESULTS_DIR`, and `COND` before running.

**Input:** `{condition}_combined.tsv`
**Output:** `{condition}_expressed.tsv`
**R packages:** `dplyr`, `readr`, `Segmentor3IsBack`

---

### Step 3 — UTR table and merge (R language)

**Script:** `scripts/ReferenceConstruction/03_utr_table.R`

Converts expressed-region tables to per-gene UTR lengths, takes the per-boundary maximum across conditions and replicates to produce `DRS_UTR_corrected.tsv`, then merges with `Nagalakshmi_UTR.csv` to produce `final_utr.tsv`. A built-in sanity check confirms zero merge violations (i.e. no final boundary is shorter than either input).

**Input:** `{cond}_expressed.tsv` for all conditions; `Nagalakshmi_UTR.csv`
**Output:** `DRS_UTR_corrected.tsv`, `final_utr.tsv` (and intermediate per-condition / per-group files)
**R packages:** `dplyr`, `readr`

---

### Step 4 — Build transcript FASTA reference (Python)

**Script:** `scripts/ReferenceConstruction/04_build_utr_reference.py`

Reconstructs a spliced mRNA FASTA from the UTR-length table, the curated gene-model backbone GFF3, and the SGD genomic FASTA. Internally: (1) renames NCBI chromosome headers to bare roman numerals; (2) overlays UTR coordinates onto the GFF3; (3) calls `gffread` to extract spliced transcripts.

```bash
python scripts/ReferenceConstruction/04_build_utr_reference.py \
    --template-gff  Data/Annotation/gene_models_backbone_with_UTR_slots.gff3 \
    --utr-table     Data/Annotation/final_utr.tsv \
    --genome        Data/Annotation/S288C_reference_sequence_R64-4-1_20230823.fsa \
    --out-prefix    FinalReference/final
```

**Dependencies:** Python ≥ 3.8, `gffread` (in `$PATH`)

---

### Preprocessing — alignment and metrics (HPC shell scripts)

These scripts are written for PBS/Torque HPC systems (tested on NCI Gadi). All paths and scheduler directives are parameterised in a `CONFIGURATION` block at the top of each script — edit these before submitting. Scripts 02–05 depend on output from earlier scripts.

| Script | Description |
|---|---|
| `01_build_STAR_indices.sh` | Builds STAR indices for the four transcriptome references (ORF-only, ORF±1000, our_ref, Nagalakshmi). Run before `03_align_ShortRead.sh`. |
| `02_align_DRS.sh` | Aligns DRS reads to all five references using minimap2 (`-ax splice -uf -k14` for the genome; `-ax map-ont -k14 -N10` for transcriptomes). Outputs sorted, indexed BAMs. |
| `03_align_ShortRead.sh` | Aligns short-read data to the genome with bwa-mem2 and to the four transcriptome references with STAR. Merges per-sample BAMs with `samtools merge`. |
| `04_metagene_coverage.sh` | deepTools-based metagene pipeline (BAM → bigWig → `computeMatrix`). Retained for reproducibility; the pysam-based approach in `05_metrics.sh` Section 5 is used for the paper figures. |
| `05_metrics.sh` | All post-alignment metrics in one job: (1) true SR mapping rate from FASTQ read counts; (2) flagstat, multimapping, antisense, and alignments-per-read for all BAMs; (3) MAPQ=0 multi-mapping proxy for bwa-mem2 genome BAMs; (4) soft-clip statistics for DRS BAMs; (5) pysam-based metagene profiles restricted to the shared gene set. |
| `compute_metagene.py` | Called by `05_metrics.sh` Section 5. Computes mean ± SEM coverage profiles anchored at the ORF start and end for any transcriptome BAM. Accepts per-gene UTR lengths or fixed flanks, and can restrict all references to the same shared gene set for a fair comparison. |
| `softclip_metrics.py` | Called by `05_metrics.sh` Section 4. Summarises soft-clip length distributions at the 5′ and 3′ ends of DRS alignments. |

> **Note on multi-mapping:** bwa-mem2 never writes secondary alignment records (SAM flag `0x100`) regardless of alignment flags. Multi-mapping rates for genome BAMs are therefore computed as the fraction of primary reads with MAPQ = 0 (`05_metrics.sh` Section 3), not via secondary record counts.

> **Note on TSS/TES:** internal code, comments and file names (`compute_metagene.py`, `04_metagene_coverage.sh`, `*_TSS.tsv` / `*_TES.tsv`) use TSS and TES as shorthand for the anchor points these profiles actually use: the annotated ORF start and ORF end, not the transcript start or end. The manuscript and figures say "ORF start"/"ORF end" throughout, per Reviewer 2's request to hold that distinction; the shorthand survives here only in code-level naming.

---

### Figures

Figure plotting scripts are in `scripts/Figures/`. They read the pre-computed metric files from `Data/Metrics/` and the annotation tables from `Data/Annotation/`, so they can be run locally without re-running the HPC pipeline.

Supplementary figure numbers below are the manuscript's current numbering; the filenames still carry an earlier numbering from before a later reordering and were not renamed to match, so don't infer the figure number from the filename.

| Script | Figure |
|---|---|
| `plot_figure2.py` | DRS UTR length distributions and comparison with Nagalakshmi *et al.* |
| `plot_figure3.py` | Alignment quality metrics across five references × two data types. |
| `plot_figure4.py` | Metagene coverage profiles anchored at the ORF start and end for all references (the axes and this table both say ORF start/end, not TSS/TES — see the note under `compute_metagene.py`). |
| `plot_suppl_fig3.R` | **Supplementary Fig. S5.** Four-panel genome-wide comparison of the pre-merge DRS calls against Nagalakshmi *et al.* |
| `plot_suppl_fig4.py` | **Supplementary Fig. S8.** Read pileups for DRS and short-read BAMs together, one gene per panel group, validating boundaries against two independent datasets neither used in construction. |
| `plot_suppl_fig2.py` | A read-pileup design for the eight-loci figure, superseded by `plot_figure1b_and_suppl_fig3.py` below once Reviewer 1 asked for genomic coordinates and a CDS track instead; kept for reference, not the current source of any released figure. |

**Figure 1 and Supplementary Fig. S3** — the segmentation worked example, redrawn on genomic coordinates with a CDS track and the Nagalakshmi annotation overlaid (Reviewer 1's specific request). Figure 1B follows one locus (YBL091C) through four steps; Supplementary Fig. S3 repeats the coverage + annotation pair for eight further loci chosen to be difficult, not representative.

| Script | Description |
|---|---|
| `locus_windows.py` | Defines each window (ORF span ± 1000 nt) and pulls every CDS/UTR/intron feature from both GFF3s that overlaps it, from `Final_UTRs.gff3` and `Nagalakshmi_UTRs.gff3`. |
| `locus_coverage_real.py` | Real per-position coverage for each window, summed over the four fraction-level conditions; the antisense track is a neighbouring gene's own sense signal, reconstructed from its own coverage table and mapped into this window's coordinates. |
| `locus_segmentation.R` | Segments each window's coverage with `segmentation_corrected.R` (see Sensitivity below), on the gene's own transcript-direction ordering rather than genomic order — the two run in opposite directions for a minus-strand gene. |
| `plot_figure1b_and_suppl_fig3.py` | Draws Figure 1B and Supplementary Fig. S3 from the three tables above. |
| `plot_figure1a_schematic.py` | Figure 1 panel A: a synthetic schematic of the segmentation pipeline. Nothing in it is measured; the caption says so. |
| `compose_figure1.py` | Stacks panel A above panel B into the final Figure 1, by translation only (no scaling, no rasterising). |

```bash
python3 scripts/Figures/locus_windows.py --annotation-dir Data/Annotation --outdir Data/LocusTracks
python3 scripts/Figures/locus_coverage_real.py --indir Data/LocusTracks \
    --combined-dir Data/CombinedProfiles --outdir Data/LocusTracks
Rscript scripts/Figures/locus_segmentation.R --indir Data/LocusTracks --outdir Data/LocusTracks
python3 scripts/Figures/plot_figure1b_and_suppl_fig3.py --indir Data/LocusTracks --outdir Figures_out
python3 scripts/Figures/plot_figure1a_schematic.py --outdir Figures_out
python3 scripts/Figures/compose_figure1.py --panel-a Figures_out/figure1a.pdf \
    --panel-b Figures_out/figure1b.pdf --out Figures_out/figure1.pdf
```

**Supplementary Fig. S9** — metagene coverage split by whether this study moved the boundary: the same our_ref genes and profiles as Figure 4, split on whether the final annotation extended that boundary by more than 20 nt.

| Script | Description |
|---|---|
| `compute_metagene_scaled.py` | Per-gene-scaled metagene profiles from the our_ref BAM, split into extended/unchanged groups, with structural zeros kept as NaN rather than padded (a mean of per-gene ratios is otherwise dominated by the smallest denominators). |
| `plot_suppl_fig9.py` | Draws Supplementary Fig. S9 from the profiles above. |

```bash
python3 scripts/Figures/compute_metagene_scaled.py --bam DRS_our_ref.bam --data-type DRS \
    --final-utr Data/Annotation/final_utr.tsv --nagalakshmi-utr Data/Annotation/Nagalakshmi_UTR.csv \
    --restrict-genes Data/Metrics/profiles/shared_genes.txt --outdir Data/Metrics/profiles_scaled
# repeat with --bam SR_merged_our_ref.bam --data-type SR --append
python3 scripts/Figures/plot_suppl_fig9.py --profiles Data/Metrics/profiles_scaled/metagene_scaled_profiles.tsv --outdir Figures_out
```

Supplementary Fig. S4 (loci where the pre-merge DRS call is shorter than the prior reference) has no plotting script in this repository yet.

---

### Comparisons — validation against independent annotations and references

Scripts in `scripts/Comparisons/` benchmark the final annotation against data and annotations that took no part in building it.

**TIF-seq (Pelechano *et al.*, 2013)**

| Script | Description |
|---|---|
| `derive_pelechano_utr.py` | Derives `Data/Annotation/Pelechano_UTR.csv` from the raw TIF-seq isoform table (GEO GSE39128) and the gene-model backbone GFF3, selecting the major covering isoform per gene. |
| `plot_final_vs_tifseq_comparison.R` | Compares `final_utr.tsv` against `Pelechano_UTR.csv` with signed differences in both directions (the final annotation can be longer *or* shorter than TIF-seq, unlike the Nagalakshmi comparison it was built from). Reports Spearman correlation and over-/under-extension rates at a configurable threshold. |

```bash
Rscript scripts/Comparisons/plot_final_vs_tifseq_comparison.R \
    --final_utr Data/Annotation/final_utr.tsv \
    --ref_utr   Data/Annotation/Pelechano_UTR.csv \
    --outdir    Figures_out \
    --threshold 20
```

**TIF-seq, both directions explicit** — Reviewer 2 asked for the TIF-seq validation reported with explicit over-extension *and* under-extension rates, not one tail folded into "otherwise concordant or longer."

| Script | Description |
|---|---|
| `tifseq_rates.py` | Over-/under-extension at three tolerances for three annotations (final, Nagalakshmi, DRS-only); the tolerance is derived from the benchmark's own dispersion rather than asserted; a provenance breakdown and the envelope test (whether a boundary falls inside the range of ends TIF-seq actually observed). |
| `plot_suppl_fig6.R` | **Supplementary Fig. S6.** Both tails at 20 nt for three annotations, both rates as a function of tolerance, the benchmark's own dispersion, and the envelope test. |

```bash
python3 scripts/Comparisons/tifseq_rates.py --annotation-dir Data/Annotation \
    --pelechano-raw GSE39128_tsedall.txt.gz --outdir Data/TIFseq
Rscript scripts/Comparisons/plot_suppl_fig6.R --indir Data/TIFseq --outdir Figures_out
```

**Is the DRS-shorter shortfall real, or a coverage artefact?** For genes where the DRS call is shorter than the prior reference, TIF-seq arbitrates directly, and a mechanistic prediction follows from DRS being read 3′ to 5′: a processivity artefact must grow with transcript length and shrink with depth at the 5′ end, and do neither at the 3′ end.

| Script | Description |
|---|---|
| `shortfall.py` | Where DRS calls shorter, which call TIF-seq agrees with (two-way and three-way), logistic models of the shortfall against CDS length and depth, and what the merge rule's cost looks like against TIF-seq. |
| `read_geometry.py` | The same prediction tested directly on the alignments, no annotation quantity involved: how far each read's aligned start/end sits from the annotated boundary, and whether the fraction of reads reaching an end decays with transcript length. |
| `plot_suppl_fig7.R` | **Supplementary Fig. S7.** Three-way TIF-seq arbitration, the reversed-direction control, the odds-ratio mechanism from the annotation, and the same contrast read directly off the alignments. |

```bash
python3 scripts/Comparisons/shortfall.py --annotation-dir Data/Annotation \
    --softclip-dir Data/Metrics/softclip --outdir Data/Shortfall
python3 scripts/Comparisons/read_geometry.py --bam DRS_our_ref_validation.bam \
    --shortfall-genes Data/Shortfall/shortfall_genes.tsv --outdir Data/Shortfall
Rscript scripts/Comparisons/plot_suppl_fig7.R --indir Data/Shortfall --outdir Figures_out
```

`read_geometry.py` needs the validation library's own BAM (independent of the six construction libraries, so agreement with a discarded DRS call is replication, not circularity) — not the our_ref BAM used elsewhere in this repository for the pooled construction data. Running it against the wrong BAM will still produce output but the headline percentages will not match the manuscript.

**UNAGI (Al Kadi *et al.*, 2020)** — an independent nanopore full-length cDNA annotation: different library chemistry, basecaller, assembler and laboratory.

| Script | Description |
|---|---|
| `unagi_compare.py` | Gene-level comparison against UNAGI's two supplementary tables, matched to the backbone by ORF coordinate (not by accession map). Reports coverage, over-/under-extension rates at three tolerances, whether the merge rule improved or only moved agreement, a provenance breakdown, UNAGI's own definition and dispersion, and UNAGI vs. TIF-seq as two independent benchmarks measured against each other. |

The two source files (MOESM5, MOESM6) are the publisher-hosted supplementary tables of the UNAGI paper and are not redistributed here; the module docstring gives the two-line `curl` command to fetch them.

```bash
python3 scripts/Comparisons/unagi_compare.py \
    --annotation-dir Data/Annotation \
    --unagi-dir      Data/UNAGI \
    --outdir         Data/UNAGI
```

**Newly overlapping neighbours** — extending UTR boundaries in a compact genome creates new overlaps between neighbouring genes, which matters for read assignment when the overlap reaches into a neighbour's coding sequence.

| Script | Description |
|---|---|
| `neighbour_overlap.py` | Compares `Final_UTRs.gff3` against `Nagalakshmi_UTRs.gff3` directly (the two GFF3s gffread built each transcriptome FASTA from) to find gene pairs that overlap under the released annotation but did not under the prior one, and flags the subset reaching into a neighbour's ORF body on the same strand — the "worst tier" `abundance_concordance.py` tests for a measurable effect. Does not classify genes by SGD ORF status or check overlap against non-mRNA features (tRNA, snoRNA, etc.); that is a separate analysis with its own external SGD annotation input. |

```bash
python3 scripts/Comparisons/neighbour_overlap.py \
    --annotation-dir Data/Annotation --outdir Data/Overlap
```

**Abundance and read migration** — does changing the reference move per-gene abundance estimates, and does it move reads across newly overlapping gene pairs specifically?

| Script | Description |
|---|---|
| `abundance_concordance.py` | Best-hit read counts for two technologies (direct RNA, short read) under both references, compared as CPM and TPM: Spearman/Pearson correlation, fraction within 1.2-fold, fraction moving 2-fold or more, and whether the newly-overlapping "worst tier" genes (from `neighbour_overlap.py`) move more than the rest (Cliff's delta). |
| `read_migration.sh` | Aligns the same reads to both references in one run, records each read's best-hit transcript under each, and joins the two into a transition table (`old_transcript`, `new_transcript`, `n_reads`). |
| `read_migration_analyse.py` | Classifies every transition (unchanged / gained / lost / migrated) and splits migrated reads into those attributable to a newly overlapping pair versus everything else, separating out the RDN37 rDNA repeat, a known multi-mapping sink independent of this annotation. |

```bash
scripts/Comparisons/read_migration.sh \
    Data/Annotation/Final_w.fa Data/Annotation/Nagalakshmi_w.fa \
    drs.fq.gz sr.fq.gz work/migration

python3 scripts/Comparisons/read_migration_analyse.py \
    --transitions-dir work/migration --overlap-dir Data/Overlap --outdir Data/Abundance
```

`abundance_concordance.py` needs each technology's best-hit counts and a transcript-length table per reference; see its module docstring for the `minimap2`/`samtools`/`seqkit` one-liners.

---

## Sensitivity — segmentation threshold sweep

Scripts in `scripts/Sensitivity/` answer why the four segmentation thresholds (coverage floor, breakpoint-merge window, breakpoint effect size, expressed-segment floor) take the values they do, by re-running boundary calling across a 5-value-per-threshold grid (625 combinations) on the six real construction conditions.

| Script | Description |
|---|---|
| `segmentation_corrected.R` | The boundary-calling algorithm (`select_segments`, `detect_expressed_region_corrected`, `compute_utrs`) factored out of `ReferenceConstruction/02_transcripts_segmentation.R` so it can be re-run with non-released threshold values. Sourced by `sensitivity_sweep.R`, not reimplemented elsewhere. |
| `sensitivity_sweep.R` | Runs the 625-combination grid on each of the six conditions, reusing one Segmentor3IsBack change-point pass per gene per condition across the whole grid. The expensive step (~1 hour for all six conditions). |
| `sensitivity_analyse.py` | Compares every combination's six-condition-max call against the released combination: a whole-grid summary (how far the worst-case joint move across all four thresholds can shift the calls) and a per-parameter one-step breakdown (which threshold drives the movement, and how much of the gene population). |
| `sensitivity_curves.py` | **Supplementary Fig. S1.** One line per threshold, the boundary-change rate across all five grid points as that threshold moves with the other three held at the released value — the shape of the response, not just the first step. |

```bash
Rscript scripts/Sensitivity/sensitivity_sweep.R
python3 scripts/Sensitivity/sensitivity_analyse.py
python3 scripts/Sensitivity/sensitivity_curves.py
```

`sensitivity_sweep.R` reads the six `{condition}_combined.tsv` files `ReferenceConstruction/01_aggregate_data.R` produces.

---

## Dependencies

### R (reference construction)
- R ≥ 4.1
- `dplyr`, `readr` (CRAN)
- `Segmentor3IsBack` (CRAN)

### Python (reference construction, metagene, figures)
- Python ≥ 3.8
- `numpy`, `pandas`, `pysam`, `matplotlib`, `seaborn`

### Command-line tools (HPC preprocessing)
- `minimap2` ≥ 2.24
- `STAR` ≥ 2.7
- `bwa-mem2` ≥ 2.2
- `samtools` ≥ 1.16
- `gffread` ≥ 0.12 (for Step 4)
- `deepTools` ≥ 3.5 (optional; `04_metagene_coverage.sh` only)
- `mosdepth` ≥ 0.3 (optional; `05_metrics.sh`)

---

## Quick start

To reproduce the final transcript FASTA from the provided annotation files (no HPC required):

```bash
# 1. Clone the repository
git clone https://github.com/Arnaroo/Yeast_Annotation.git
cd Yeast_Annotation

# 2. Install Python dependencies
pip install numpy pandas pysam matplotlib seaborn

# 3. Build the transcript FASTA reference
python scripts/ReferenceConstruction/04_build_utr_reference.py \
    --template-gff  Data/Annotation/gene_models_backbone_with_UTR_slots.gff3 \
    --utr-table     Data/Annotation/final_utr.tsv \
    --genome        Data/Annotation/S288C_reference_sequence_R64-4-1_20230823.fsa \
    --out-prefix    FinalReference/final
# Output: FinalReference/final_w.fa  (transcript FASTA)
#         FinalReference/final_UTRs.gff3
#         FinalReference/final_genome_chr_id.fa

# 4. Plot figures (reads pre-computed metrics from Data/Metrics/)
python scripts/Figures/plot_figure2.py
python scripts/Figures/plot_figure3.py
python scripts/Figures/plot_figure4.py
```

To reproduce the annotation from raw DRS coverage profiles, run scripts 01–03 in `scripts/ReferenceConstruction/` in order, editing the `CONFIGURATION` block in each script to point to your data. The HPC alignment and metric scripts in `scripts/Preprocessing/` can then be used to reproduce the alignment benchmarking results.

---

## Citation

If you use this annotation or code, please cite:

> Rossini O, Cleynen A, Shirokikh N. A refined *Saccharomyces cerevisiae* reference transcriptome from Direct RNA Sequencing, with a reusable pipeline for UTR annotation updates. *FEMS Yeast Research*, in revision.

and, for the data and code archive itself:

> Rossini O, Cleynen A, Shirokikh N. (2026). A refined *Saccharomyces cerevisiae* reference transcriptome from Direct RNA Sequencing [data set]. Zenodo. https://doi.org/10.5281/zenodo.20828039

The Nagalakshmi *et al.* reference used in the merge:

> Nagalakshmi U, Wang Z, Waern K, Shou C, Raha D, Gerstein M, Snyder M (2008). The transcriptional landscape of the yeast genome defined by RNA sequencing. *Science* 320:1344–1349. https://doi.org/10.1126/science.1158441

---

## License

Code is released under the [MIT License](LICENSE). Annotation data files are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
