#!/bin/bash
# ============================================================
# 02_align_DRS.sh
#
# Align Direct RNA Sequencing (DRS / Oxford Nanopore) reads
# to five references using minimap2:
#   - Genome (splice-aware, -ax splice -uf)
#   - ORF only (transcriptome, -ax map-ont)
#   - ORF ±1000 nt (transcriptome, -ax map-ont)
#   - Our merged annotation (transcriptome, -ax map-ont)
#   - Nagalakshmi et al. annotation (transcriptome, -ax map-ont)
#
# Output: one sorted, indexed BAM per reference in OUTPUT_DIR.
#
# Dependencies: minimap2, samtools
# ============================================================

# ── PBS / SLURM header ────────────────────────────────────────────────────
#PBS -P YOUR_PROJECT_CODE
#PBS -q normal
#PBS -N align_DRS
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l walltime=10:00:00
#PBS -o /path/to/logs/align_DRS.out
#PBS -e /path/to/logs/align_DRS.err
##PBS -l storage=gdata/YOUR_PROJECT+scratch/YOUR_PROJECT

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these paths before submitting
# ════════════════════════════════════════════════════════════════════════

# Input DRS reads (FASTQ, can be gzipped)
DRS_READS=/path/to/DRS_reads.fastq.gz

# Reference FASTA files
REF_DIR=/path/to/references
GENOME=${REF_DIR}/genome.fa
ORF_ONLY=${REF_DIR}/orf_genomic_all.fasta
ORF_1000=${REF_DIR}/orf_genomic_1000_all.fasta
OUR_REF=${REF_DIR}/our_reference.fa
NAGALAKSHMI=${REF_DIR}/nagalakshmi_reference.fa

# Output directory for BAM files
OUTPUT_DIR=/path/to/alignments/DRS

# Number of threads (match PBS ncpus above)
THREADS=48

# ════════════════════════════════════════════════════════════════════════

module load samtools
module load minimap2

mkdir -p ${OUTPUT_DIR}

# ── Helper: align, sort, index ────────────────────────────────────────────
align_and_sort() {
    local MODE=$1       # minimap2 preset (e.g. splice, map-ont)
    local EXTRA=$2      # additional minimap2 flags (or "")
    local REF=$3        # reference FASTA
    local OUT=$4        # output BAM path

    echo "$(date) -- Aligning to $(basename ${REF} .fa)$(basename ${REF} .fasta)..."
    minimap2 -ax ${MODE} ${EXTRA} \
        -t ${THREADS} \
        ${REF} ${DRS_READS} \
    | samtools sort -@ ${THREADS} -o ${OUT}
    samtools index ${OUT}
    echo "$(date) -- Done: $(basename ${OUT})"
}

# ── Genome (splice-aware, strand-sensitive for DRS poly-A) ────────────────
# -uf : force strand-specific alignment (U = poly-A captured, + strand)
# -k14: shorter k-mer for noisy DRS reads
align_and_sort splice "-uf -k14" \
    ${GENOME} ${OUTPUT_DIR}/DRS_genome.bam

# ── Transcriptome references (map-ont mode) ───────────────────────────────
# -k14 : shorter k-mer for noisy DRS reads
# -N10 : retain up to 10 secondary alignments (for multi-mapping analysis)

align_and_sort map-ont "-k14 -N10" \
    ${ORF_ONLY}    ${OUTPUT_DIR}/DRS_orf_only.bam

align_and_sort map-ont "-k14 -N10" \
    ${ORF_1000}    ${OUTPUT_DIR}/DRS_orf_1000.bam

align_and_sort map-ont "-k14 -N10" \
    ${OUR_REF}     ${OUTPUT_DIR}/DRS_our_ref.bam

align_and_sort map-ont "-k14 -N10" \
    ${NAGALAKSHMI} ${OUTPUT_DIR}/DRS_nagalakshmi.bam

echo "$(date) -- All DRS alignments complete."
echo "$(date) -- Output directory: ${OUTPUT_DIR}"
