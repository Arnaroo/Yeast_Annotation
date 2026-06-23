#!/bin/bash
# ============================================================
# 01_build_STAR_indices.sh
#
# Build STAR indices for all transcriptome references used in
# short-read alignment (script 03_align_ShortRead.sh).
#
# STAR is used for transcriptome references only. For the
# genome reference, bwa-mem2 is used (no GTF required).
#
# For transcriptome FASTA references, --genomeSAindexNbases is
# set to 11 (formula: min(14, log2(GenomeLength)/2 - 1) for
# ~14 Mb yeast transcriptomes). Adjust if STAR warns.
#
# Run this script BEFORE 03_align_ShortRead.sh.
# ============================================================

# ── PBS / SLURM header ────────────────────────────────────────────────────
# Adjust resource requests and scheduler directives for your HPC system.
#PBS -P YOUR_PROJECT_CODE
#PBS -q normal
#PBS -N build_STAR_indices
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l walltime=04:00:00
# Edit log paths and storage flags for your system:
#PBS -o /path/to/logs/build_STAR_indices.out
#PBS -e /path/to/logs/build_STAR_indices.err
##PBS -l storage=gdata/YOUR_PROJECT+scratch/YOUR_PROJECT

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these paths before submitting
# ════════════════════════════════════════════════════════════════════════

# Path to the STAR executable (or use $(which STAR) if in PATH)
STAR=$(which STAR)

# Directory containing reference FASTA files
REF_DIR=/path/to/references

# Directory where STAR indices will be written
INDEX_DIR=${REF_DIR}/STAR_indices

# Number of threads (match PBS ncpus above)
THREADS=48

# Reference FASTA filenames — edit to match your files
FASTA_ORF_ONLY=${REF_DIR}/orf_genomic_all.fasta
FASTA_ORF_1000=${REF_DIR}/orf_genomic_1000_all.fasta
FASTA_OUR_REF=${REF_DIR}/our_reference.fa
FASTA_NAGALAKSHMI=${REF_DIR}/nagalakshmi_reference.fa

# ════════════════════════════════════════════════════════════════════════

module load samtools

mkdir -p ${INDEX_DIR}

# ── ORF only ──────────────────────────────────────────────────────────────
echo "$(date) -- Building STAR index: ORF_ONLY..."
mkdir -p ${INDEX_DIR}/orf_only
${STAR} --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --genomeDir ${INDEX_DIR}/orf_only \
    --genomeFastaFiles ${FASTA_ORF_ONLY} \
    --genomeSAindexNbases 11 \
    --genomeChrBinNbits 11
echo "$(date) -- Done: orf_only"

# ── ORF ±1000 nt ──────────────────────────────────────────────────────────
echo "$(date) -- Building STAR index: ORF_1000..."
mkdir -p ${INDEX_DIR}/orf_1000
${STAR} --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --genomeDir ${INDEX_DIR}/orf_1000 \
    --genomeFastaFiles ${FASTA_ORF_1000} \
    --genomeSAindexNbases 11 \
    --genomeChrBinNbits 11
echo "$(date) -- Done: orf_1000"

# ── Our merged annotation ─────────────────────────────────────────────────
echo "$(date) -- Building STAR index: OUR_REF..."
mkdir -p ${INDEX_DIR}/our_ref
${STAR} --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --genomeDir ${INDEX_DIR}/our_ref \
    --genomeFastaFiles ${FASTA_OUR_REF} \
    --genomeSAindexNbases 11 \
    --genomeChrBinNbits 11
echo "$(date) -- Done: our_ref"

# ── Nagalakshmi et al. annotation ────────────────────────────────────────
echo "$(date) -- Building STAR index: NAGALAKSHMI..."
mkdir -p ${INDEX_DIR}/nagalakshmi
${STAR} --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --genomeDir ${INDEX_DIR}/nagalakshmi \
    --genomeFastaFiles ${FASTA_NAGALAKSHMI} \
    --genomeSAindexNbases 11 \
    --genomeChrBinNbits 11
echo "$(date) -- Done: nagalakshmi"

echo "$(date) -- All STAR indices built successfully."
echo "$(date) -- Index directory: ${INDEX_DIR}"
