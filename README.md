# An updated UTR annotation for *Saccharomyces cerevisiae* from Direct RNA Sequencing

This repository contains the data, code, and pre-built reference files associated with:

> **Rossini _et al._** (2026). An updated UTR annotation for *Saccharomyces cerevisiae* derived from Oxford Nanopore Direct RNA Sequencing.

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
│   ├── Bam/                 # Example BAM files (DRS and short-read)
│   └── Metrics/             # Pre-computed alignment quality metrics
│       ├── aln_per_read/
│       ├── antisense/
│       ├── flagstat/
│       ├── profiles/        # Metagene coverage profiles (TSV)
│       └── softclip/
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
│   └── Figures/                 # Figure plotting scripts
│       ├── plot_figure2.py
│       ├── plot_figure3.py
│       ├── plot_figure4.py
│       ├── plot_suppl_fig2.py
│       ├── plot_suppl_fig3.R
│       └── plot_suppl_fig4.py
└── README.md
```

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

### BAM files (`Data/Bam/`)

Representative sorted, indexed BAM files illustrating alignments to the final reference:

| File | Description |
|---|---|
| `DRS_our_ref.bam` | DRS reads (minimap2, map-ont) aligned to `Final_w.fa`. |
| `SR_merged_our_ref.bam` | Short-read data (STAR) aligned to `Final_w.fa`. |
| `NS_PS_Rep3_filtered.sorted.bam` | Filtered DRS reads for the NS_PS replicate 3 condition (used in annotation construction). |

### Metrics (`Data/Metrics/`)

Pre-computed alignment quality metrics for all five references × two data types, ready for direct use with the figure scripts:

- `flagstat/` — `samtools flagstat` output and multimapping rates
- `antisense/` — antisense alignment rates
- `aln_per_read/` — alignments-per-read distributions
- `softclip/` — soft-clip statistics (DRS only)
- `profiles/` — metagene coverage profiles around the ORF start (TSS) and end (TES)

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
| `compute_metagene.py` | Called by `05_metrics.sh` Section 5. Computes mean ± SEM coverage profiles anchored at the true ORF start (TSS) and end (TES) for any transcriptome BAM. Accepts per-gene UTR lengths or fixed flanks, and can restrict all references to the same shared gene set for a fair comparison. |
| `softclip_metrics.py` | Called by `05_metrics.sh` Section 4. Summarises soft-clip length distributions at the 5′ and 3′ ends of DRS alignments. |

> **Note on multi-mapping:** bwa-mem2 never writes secondary alignment records (SAM flag `0x100`) regardless of alignment flags. Multi-mapping rates for genome BAMs are therefore computed as the fraction of primary reads with MAPQ = 0 (`05_metrics.sh` Section 3), not via secondary record counts.

---

### Figures

Figure plotting scripts are in `scripts/Figures/`. They read the pre-computed metric files from `Data/Metrics/` and the annotation tables from `Data/Annotation/`, so they can be run locally without re-running the HPC pipeline.

| Script | Figure |
|---|---|
| `plot_figure2.py` | DRS UTR length distributions and comparison with Nagalakshmi *et al.* |
| `plot_figure3.py` | Alignment quality metrics across five references × two data types. |
| `plot_figure4.py` | Metagene coverage profiles at TSS and TES for all references. |
| `plot_suppl_fig2.py` | Supplementary: four-panel DRS vs. Nagalakshmi UTR comparison. |
| `plot_suppl_fig3.R` | Supplementary: per-condition UTR length distributions. |
| `plot_suppl_fig4.py` | Supplementary: soft-clip length distributions. |

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
git clone https://github.com/Arnaroo/yeast-utr-annotation.git
cd yeast-utr-annotation

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

> Rossini *et al.* (2026). An updated UTR annotation for *Saccharomyces cerevisiae* derived from Oxford Nanopore Direct RNA Sequencing. 

The Nagalakshmi *et al.* reference used in the merge:

> Nagalakshmi U, Wang Z, Waern K, Shou C, Raha D, Gerstein M, Snyder M (2008). The transcriptional landscape of the yeast genome defined by RNA sequencing. *Science* 320:1344–1349. https://doi.org/10.1126/science.1158441

---

## License

Code is released under the [MIT License](LICENSE). Annotation data files are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
