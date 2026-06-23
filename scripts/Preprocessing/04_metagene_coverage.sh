#!/bin/bash
# ============================================================
# 04_metagene_coverage.sh
#
# Compute metagene (reference-point) coverage profiles around
# ORF start (TSS) and ORF end (TES) boundaries for all five
# references, using the deepTools bamCoverage + computeMatrix
# pipeline.
#
# Workflow per reference:
#   1. BAM → bigWig (CPM-normalised, bamCoverage)
#   2. BED from FASTA index (transcript coordinates)
#   3. computeMatrix reference-point at TSS and TES (±1000 bp)
#
# All references use the same transcript universe for fair
# comparison. Each BAM is paired only with its own BED to avoid
# chromosome namespace mismatches (orf_only/orf_1000 use bare
# gene names; our_ref/nagalakshmi use gene_mRNA suffix).
#
# Key flags:
#   --missingDataAsZero : fill empty bins with 0 (not NaN)
#                         avoids empty-matrix crash for orf_only
#   bamCoverage --samFlagExclude 2304 : skip secondary (0x100)
#                                       and supplementary (0x800)
#
# Note: this deepTools pipeline was developed and partially
# superseded by the pysam-based compute_metagene.py approach
# (run_metagene_profiles.sh), which is used for Figure 4.
# This script is retained for reproducibility.
#
# Dependencies: samtools, bamCoverage (deepTools), computeMatrix
#               (deepTools), pyBigWig (Python)
# ============================================================

# ── PBS / SLURM header ────────────────────────────────────────────────────
#PBS -P YOUR_PROJECT_CODE
#PBS -q normal
#PBS -N metagene_coverage
#PBS -l mem=128GB
#PBS -l ncpus=16
#PBS -l wd
#PBS -l walltime=10:00:00
#PBS -o /path/to/logs/metagene_coverage.out
#PBS -e /path/to/logs/metagene_coverage.err
##PBS -l storage=gdata/YOUR_PROJECT+scratch/YOUR_PROJECT

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these paths before submitting
# ════════════════════════════════════════════════════════════════════════

# DRS and short-read BAM directories (output of 02 and 03)
DRS_DIR=/path/to/alignments/DRS
SR_DIR=/path/to/alignments/ShortRead

# Reference FASTA files (used to build BED files)
REF_DIR=/path/to/references
FASTA_ORF_ONLY=${REF_DIR}/orf_genomic_all.fasta
FASTA_ORF_1000=${REF_DIR}/orf_genomic_1000_all.fasta
FASTA_OUR_REF=${REF_DIR}/our_reference.fa
FASTA_NAGALAKSHMI=${REF_DIR}/nagalakshmi_reference.fa

# Output directories
OUT_DIR=/path/to/metagene_output
BW_DIR=${OUT_DIR}/bigwig
BED_DIR=${OUT_DIR}/beds
MAT_DIR=${OUT_DIR}/matrices

# Total primary SR read count (used for CPM normalisation).
# Obtain with: samtools view -F 0x904 -c SR_merged_genome.bam
SR_TOTAL=FILL_IN_YOUR_VALUE

# Number of threads (match PBS ncpus above)
THREADS=16

# ════════════════════════════════════════════════════════════════════════

module load python3
module load samtools
export PATH=$HOME/.local/bin:$PATH   # ensure deepTools binaries are found

mkdir -p ${BW_DIR} ${BED_DIR} ${MAT_DIR}

set -euo pipefail

# ── Library size scale factors (CPM) ─────────────────────────────────────
DRS_TOTAL=$(samtools view -@ ${THREADS} -F 0x904 -c ${DRS_DIR}/DRS_genome.bam)
echo "$(date) -- DRS total primary reads: ${DRS_TOTAL}"

SR_SCALE=$(awk  "BEGIN{printf \"%.8f\", 1000000/${SR_TOTAL}}")
DRS_SCALE=$(awk "BEGIN{printf \"%.8f\", 1000000/${DRS_TOTAL}}")
echo "$(date) -- SR  CPM scale factor: ${SR_SCALE}"
echo "$(date) -- DRS CPM scale factor: ${DRS_SCALE}"

# ── STEP 1: Generate BED files from FASTA indices ─────────────────────────
# Each BED represents all transcripts in one reference.
# BED format: chrom  start  end  name  score  strand
# Strand is set to "+" (all transcripts are sense-strand in these
# transcriptome FASTA files; strand is not used by computeMatrix
# reference-point mode).

make_bed_from_fai() {
    local FASTA=$1 BED=$2
    [ -f ${FASTA}.fai ] || samtools faidx ${FASTA}
    if [ ! -f ${BED} ]; then
        echo "$(date) -- Generating BED: $(basename ${BED})"
        awk 'BEGIN{OFS="\t"}{print $1, 0, $2, $1, ".", "+"}' \
            ${FASTA}.fai > ${BED}
    else
        echo "$(date) -- BED exists: $(basename ${BED})"
    fi
}

make_bed_from_fai ${FASTA_ORF_ONLY}    ${BED_DIR}/orf_only.bed
make_bed_from_fai ${FASTA_ORF_1000}    ${BED_DIR}/orf_1000.bed
make_bed_from_fai ${FASTA_OUR_REF}     ${BED_DIR}/our_ref.bed
make_bed_from_fai ${FASTA_NAGALAKSHMI} ${BED_DIR}/nagalakshmi.bed

# ── STEP 2: BAM → bigWig (CPM-normalised) ────────────────────────────────
# --samFlagExclude 2304 : skip secondary (0x100) + supplementary (0x800)
# --scaleFactor         : CPM normalisation (see above)
# --binSize 10          : 10 bp resolution

make_bigwig() {
    local BAM=$1 BW=$2 SCALE=$3
    if [ -f ${BW} ]; then
        echo "$(date) -- bigWig exists: $(basename ${BW})"
        return
    fi
    echo "$(date) -- bamCoverage: $(basename ${BW})"
    bamCoverage \
        --bam            ${BAM} \
        --outFileName    ${BW} \
        --outFileFormat  bigwig \
        --binSize        10 \
        --numberOfProcessors ${THREADS} \
        --normalizeUsing None \
        --scaleFactor    ${SCALE} \
        --samFlagExclude 2304 \
        --ignoreDuplicates
}

# Validate bigWig has non-zero content (catches namespace mismatches early)
validate_bw() {
    local BW=$1
    python3 - << PYEOF
import pyBigWig
bw = pyBigWig.open("${BW}")
chroms = bw.chroms()
total = sum(bw.stats(c, type="sum", exact=True) or [0]
            for c in list(chroms.keys())[:10])
bw.close()
if not total:
    raise ValueError("BigWig appears empty: ${BW}")
print(f"Validated $(basename ${BW}): sum of first 10 chroms = {total:.1f}")
PYEOF
}

# SR bigWigs
make_bigwig ${SR_DIR}/SR_merged_orf_only.bam    ${BW_DIR}/SR_orf_only.bw    ${SR_SCALE}
make_bigwig ${SR_DIR}/SR_merged_orf_1000.bam    ${BW_DIR}/SR_orf_1000.bw    ${SR_SCALE}
make_bigwig ${SR_DIR}/SR_merged_our_ref.bam     ${BW_DIR}/SR_our_ref.bw     ${SR_SCALE}
make_bigwig ${SR_DIR}/SR_merged_nagalakshmi.bam ${BW_DIR}/SR_nagalakshmi.bw ${SR_SCALE}

# DRS bigWigs
make_bigwig ${DRS_DIR}/DRS_orf_only.bam    ${BW_DIR}/DRS_orf_only.bw    ${DRS_SCALE}
make_bigwig ${DRS_DIR}/DRS_orf_1000.bam    ${BW_DIR}/DRS_orf_1000.bw    ${DRS_SCALE}
make_bigwig ${DRS_DIR}/DRS_our_ref.bam     ${BW_DIR}/DRS_our_ref.bw     ${DRS_SCALE}
make_bigwig ${DRS_DIR}/DRS_nagalakshmi.bam ${BW_DIR}/DRS_nagalakshmi.bw ${DRS_SCALE}

echo "$(date) -- Validating bigWigs..."
for BW in ${BW_DIR}/*.bw; do validate_bw ${BW}; done
echo "$(date) -- All bigWigs done."

# ── STEP 3: computeMatrix reference-point ────────────────────────────────
# Anchored at ORF start (TSS) and ORF end (TES), ±1000 bp window.
# Each BAM/bigWig is paired with its own BED (same contig namespace).
#
# --missingDataAsZero : fill empty bins with 0 rather than NaN.
#   Critical for orf_only where flanking regions have no sequence
#   and bins would otherwise be NaN, causing matrix-empty crashes.
# --outFileNameMatrix : plain-text tab file read by plot_metagene.py

run_refpoint() {
    local BW=$1 BED=$2 POINT=$3 PREFIX=$4
    local OUT_GZ=${MAT_DIR}/${PREFIX}_${POINT}.mat.gz
    local OUT_TAB=${MAT_DIR}/${PREFIX}_${POINT}.tab

    echo "$(date) -- computeMatrix ${POINT}: ${PREFIX}"
    computeMatrix reference-point \
        --scoreFileName      ${BW} \
        --regionsFileName    ${BED} \
        --samplesLabel       "${PREFIX}" \
        --referencePoint     ${POINT} \
        --upstream           1000 \
        --downstream         1000 \
        --binSize            10 \
        --missingDataAsZero \
        --numberOfProcessors ${THREADS} \
        --outFileName        ${OUT_GZ} \
        --outFileNameMatrix  ${OUT_TAB}

    NROWS=$(awk 'NF==200' ${OUT_TAB} | wc -l)
    echo "$(date) -- ${PREFIX}_${POINT}: ${NROWS} transcripts with data"
    if [ "${NROWS}" -eq 0 ]; then
        echo "WARNING: empty output for ${PREFIX}_${POINT} — check BED/bigWig namespace match"
    fi
}

# SR matrices
for REF in orf_only orf_1000 our_ref nagalakshmi; do
    run_refpoint ${BW_DIR}/SR_${REF}.bw  ${BED_DIR}/${REF}.bed  TSS SR_${REF}
    run_refpoint ${BW_DIR}/SR_${REF}.bw  ${BED_DIR}/${REF}.bed  TES SR_${REF}
done

# DRS matrices
for REF in orf_only orf_1000 our_ref nagalakshmi; do
    run_refpoint ${BW_DIR}/DRS_${REF}.bw ${BED_DIR}/${REF}.bed  TSS DRS_${REF}
    run_refpoint ${BW_DIR}/DRS_${REF}.bw ${BED_DIR}/${REF}.bed  TES DRS_${REF}
done

echo "$(date) -- All matrices complete."
echo "$(date) -- Download ${MAT_DIR}/*.tab and run plot_metagene.py locally."
