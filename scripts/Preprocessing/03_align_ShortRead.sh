#!/bin/bash
# ============================================================
# 03_align_ShortRead.sh
#
# Align short-read (Illumina) data to five references:
#
#   Genome reference   — bwa-mem2 (splice-unaware; -h 10 to
#                        report up to 10 alternative loci in
#                        the XA tag for multi-mapping analysis;
#                        note bwa-mem2 does NOT write secondary
#                        records via flag 0x100 regardless of -h)
#   Transcriptome refs — STAR (splice-suppressed in transcriptome
#                        mode; NH tag used for multi-mapping)
#
# Two short-read samples are aligned separately then merged with
# samtools merge. Edit SR1 and SR2 (and add more if needed).
#
# Run 01_build_STAR_indices.sh before this script.
#
# Dependencies: bwa-mem2, STAR, samtools
# ============================================================

# ── PBS / SLURM header ────────────────────────────────────────────────────
#PBS -P YOUR_PROJECT_CODE
#PBS -q normal
#PBS -N align_ShortRead
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l walltime=10:00:00
#PBS -o /path/to/logs/align_ShortRead.out
#PBS -e /path/to/logs/align_ShortRead.err
##PBS -l storage=gdata/YOUR_PROJECT+scratch/YOUR_PROJECT

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these paths before submitting
# ════════════════════════════════════════════════════════════════════════

# Input short-read FASTQs (add SR3, SR4 … as needed)
SR1=/path/to/sample1.fastq
SR2=/path/to/sample2.fastq

# Reference files
REF_DIR=/path/to/references
GENOME=${REF_DIR}/genome.fa
INDEX_DIR=${REF_DIR}/STAR_indices   # built by 01_build_STAR_indices.sh

# Output directories
OUTPUT_DIR=/path/to/alignments/ShortRead
TMP_DIR=${OUTPUT_DIR}/tmp

# Path to STAR (or use $(which STAR) if in PATH)
STAR=$(which STAR)

# Number of threads (match PBS ncpus above)
THREADS=48

# ════════════════════════════════════════════════════════════════════════

module load samtools
module load bwa-mem2

mkdir -p ${OUTPUT_DIR} ${TMP_DIR}

# ════════════════════════════════════════════════════════════════════════
# GENOME ALIGNMENT with bwa-mem2
#
# -h 10 : report up to 10 alternative loci per read in the XA
#         auxiliary tag. Note: bwa-mem2 does NOT write secondary
#         alignment records (flag 0x100) regardless of this flag.
#         Multi-mapping is therefore measured via MAPQ=0 primary
#         reads rather than secondary record counts (see
#         compute_metrics.sh for details).
# ════════════════════════════════════════════════════════════════════════

# Build bwa-mem2 index if not already present
if [ ! -f ${GENOME}.bwt.2bit.64 ]; then
    echo "$(date) -- Building bwa-mem2 genome index..."
    bwa-mem2 index ${GENOME}
fi

echo "$(date) -- Aligning SR1 to genome..."
bwa-mem2 mem -t ${THREADS} -h 10 ${GENOME} ${SR1} \
    | samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/SR_sample1_genome.bam
samtools index ${OUTPUT_DIR}/SR_sample1_genome.bam

echo "$(date) -- Aligning SR2 to genome..."
bwa-mem2 mem -t ${THREADS} -h 10 ${GENOME} ${SR2} \
    | samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/SR_sample2_genome.bam
samtools index ${OUTPUT_DIR}/SR_sample2_genome.bam

echo "$(date) -- Merging genome BAMs..."
samtools merge -@ ${THREADS} -f \
    ${OUTPUT_DIR}/SR_merged_genome.bam \
    ${OUTPUT_DIR}/SR_sample1_genome.bam \
    ${OUTPUT_DIR}/SR_sample2_genome.bam
samtools index ${OUTPUT_DIR}/SR_merged_genome.bam
echo "$(date) -- Done: SR_merged_genome.bam"

# ════════════════════════════════════════════════════════════════════════
# TRANSCRIPTOME ALIGNMENTS with STAR
#
# STAR is run in transcriptome mode (no GTF; FASTA reference only).
# Antisense alignments are suppressed by default in this mode
# (--outSAMstrandField is not set), so the antisense rate metric
# will be N/A for STAR BAMs.
# NH tag is written and used for multi-mapping rate calculation.
# ════════════════════════════════════════════════════════════════════════
align_star() {
    local REF_NAME=$1
    local INDEX=$2

    echo "$(date) -- STAR aligning SR1 → ${REF_NAME}..."
    mkdir -p ${TMP_DIR}/star_${REF_NAME}_s1
    ${STAR} --runThreadN ${THREADS} \
        --genomeDir ${INDEX} \
        --readFilesIn ${SR1} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM \
        --outFileNamePrefix ${TMP_DIR}/star_${REF_NAME}_s1/ \
        --outTmpDir ${TMP_DIR}/star_${REF_NAME}_s1_tmp

    echo "$(date) -- STAR aligning SR2 → ${REF_NAME}..."
    mkdir -p ${TMP_DIR}/star_${REF_NAME}_s2
    ${STAR} --runThreadN ${THREADS} \
        --genomeDir ${INDEX} \
        --readFilesIn ${SR2} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM \
        --outFileNamePrefix ${TMP_DIR}/star_${REF_NAME}_s2/ \
        --outTmpDir ${TMP_DIR}/star_${REF_NAME}_s2_tmp

    echo "$(date) -- Merging ${REF_NAME} BAMs..."
    samtools merge -@ ${THREADS} -f \
        ${OUTPUT_DIR}/SR_merged_${REF_NAME}.bam \
        ${TMP_DIR}/star_${REF_NAME}_s1/Aligned.sortedByCoord.out.bam \
        ${TMP_DIR}/star_${REF_NAME}_s2/Aligned.sortedByCoord.out.bam
    samtools index ${OUTPUT_DIR}/SR_merged_${REF_NAME}.bam
    echo "$(date) -- Done: SR_merged_${REF_NAME}.bam"
}

align_star "orf_only"    ${INDEX_DIR}/orf_only
align_star "orf_1000"    ${INDEX_DIR}/orf_1000
align_star "our_ref"     ${INDEX_DIR}/our_ref
align_star "nagalakshmi" ${INDEX_DIR}/nagalakshmi

echo "$(date) -- All short-read alignments complete."
echo "$(date) -- Output directory: ${OUTPUT_DIR}"
