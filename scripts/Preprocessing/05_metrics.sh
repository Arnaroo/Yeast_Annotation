#!/bin/bash
# ============================================================
# 05_metrics.sh
#
# All post-alignment metrics in one job, run sequentially:
#
#   SECTION 1  True SR mapping rate from FASTQ read counts
#   SECTION 2  Core alignment metrics (flagstat, multimapping,
#              antisense, alignments-per-read) for all BAMs
#   SECTION 3  MAPQ=0 multi-mapping proxy for genome BAMs
#              (bwa-mem2 never writes secondary records)
#   SECTION 4  Soft-clip statistics — DRS BAMs only
#              (calls softclip_metrics.py)
#   SECTION 5  Metagene coverage profiles — all references
#              (calls compute_metagene.py)
#
# Required Python scripts (place in SCRIPTS_DIR, see below):
#   softclip_metrics.py  — Section 4
#   compute_metagene.py  — Section 5
#
# Required Python packages (install once):
#   pip install pysam numpy pandas
#
# Required tools: samtools, mosdepth (optional, Section 2)
#
# Sections can be individually disabled by setting the
# corresponding RUN_* variable to 0.
# ============================================================

# ── PBS / SLURM header ────────────────────────────────────────────────────
# Adjust resource requests and scheduler directives for your HPC system.
# Walltime is long because sections run sequentially; split into separate
# jobs if your queue has a shorter limit.
#PBS -P YOUR_PROJECT_CODE
#PBS -q normal
#PBS -N metrics
#PBS -l mem=128GB
#PBS -l ncpus=16
#PBS -l wd
#PBS -l walltime=36:00:00
#PBS -o /path/to/logs/metrics.out
#PBS -e /path/to/logs/metrics.err
##PBS -l storage=gdata/YOUR_PROJECT+scratch/YOUR_PROJECT

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit all paths in this block before submitting
# ════════════════════════════════════════════════════════════════════════

# ── BAM directories (outputs of 02_align_DRS.sh / 03_align_ShortRead.sh)
DRS_DIR=/path/to/alignments/DRS
SR_DIR=/path/to/alignments/ShortRead

# ── Raw FASTQ files (needed for true SR mapping rate, Section 1)
SR_FASTQ_1=/path/to/sample1.fastq
SR_FASTQ_2=/path/to/sample2.fastq

# ── UTR annotation files (needed for metagene profiles, Section 5)
OUR_UTR=/path/to/final_utr.tsv          # cols: Gene  final_five_prime_utr  final_three_prime_utr
NAGALAKSHMI_UTR=/path/to/reference_UTR.csv  # cols: Gene  five_prime_utr  three_prime_utr

# ── Python scripts (softclip_metrics.py and compute_metagene.py)
SCRIPTS_DIR=/path/to/analysis/scripts

# ── Output directories
METRICS_DIR=/path/to/output/Metrics
PROFILES_DIR=${METRICS_DIR}/profiles

# ── mosdepth binary path (optional; set to "" to skip mosdepth)
MOSDEPTH=$(which mosdepth 2>/dev/null || echo "")

# ── Threads (match PBS ncpus above)
THREADS=16

# ── Toggle sections (1 = run, 0 = skip)
RUN_MAPPING_RATE=1   # Section 1: true SR mapping rate
RUN_CORE_METRICS=1   # Section 2: flagstat, multimapping, antisense, aln_per_read
RUN_BWA_MAPQ0=1      # Section 3: MAPQ=0 multi-mapping proxy for genome BAMs
RUN_SOFTCLIP=1       # Section 4: soft-clip metrics (DRS only)
RUN_METAGENE=1       # Section 5: metagene coverage profiles

# ════════════════════════════════════════════════════════════════════════

module load samtools
module load python3
set -euo pipefail

mkdir -p ${METRICS_DIR}/flagstat \
         ${METRICS_DIR}/idxstats \
         ${METRICS_DIR}/mosdepth \
         ${METRICS_DIR}/antisense \
         ${METRICS_DIR}/aln_per_read \
         ${METRICS_DIR}/softclip \
         ${PROFILES_DIR}

# ════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ════════════════════════════════════════════════════════════════════════

run_flagstat() {
    local BAM=$1 NAME=$2
    echo "$(date) -- flagstat: ${NAME}"
    samtools flagstat -@ ${THREADS} ${BAM} \
        > ${METRICS_DIR}/flagstat/${NAME}.flagstat.txt
}

# mosdepth is optional — skipped if binary not found
run_mosdepth() {
    local BAM=$1 NAME=$2
    [ -z "${MOSDEPTH}" ] && { echo "  (mosdepth not found, skipping)"; return; }
    echo "$(date) -- mosdepth: ${NAME}"
    ${MOSDEPTH} \
        --threads ${THREADS} \
        --quantize 0:1:5:10:50: \
        ${METRICS_DIR}/mosdepth/${NAME} \
        ${BAM}
}

# Multi-mapping for STAR (transcriptome) BAMs — uses NH tag
run_multimapping_STAR() {
    local BAM=$1 NAME=$2
    local OUT=${METRICS_DIR}/flagstat/${NAME}.multimapping.txt
    echo "$(date) -- multimapping (NH tag, STAR): ${NAME}"
    TOTAL=$(samtools view -@ ${THREADS} -F 0x904 -c ${BAM})
    MULTI=$(samtools view -@ ${THREADS} -F 0x904 ${BAM} \
        | awk 'BEGIN{c=0}{
            for(i=12;i<=NF;i++){
                if($i~/^NH:i:/){
                    n=substr($i,6)+0
                    if(n>=2){c++;break}
                }
            }
          }END{print c}')
    echo -e "total_primary\tmulti_mapped\tmulti_rate" > ${OUT}
    echo -e "${TOTAL}\t${MULTI}\t$(echo "scale=4; ${MULTI}/${TOTAL}" | bc)" >> ${OUT}
}

# Multi-mapping for DRS (minimap2 transcriptome) BAMs — secondary flag
run_multimapping_minimap2() {
    local BAM=$1 NAME=$2
    local OUT=${METRICS_DIR}/flagstat/${NAME}.multimapping.txt
    echo "$(date) -- multimapping (secondary flag, minimap2): ${NAME}"
    PRIMARY=$(samtools view -@ ${THREADS} -F 0x904 -c ${BAM})
    SECONDARY=$(samtools view -@ ${THREADS} -f 0x100 -F 0x804 -c ${BAM})
    echo -e "primary_mapped\tsecondary_alignments\tsecondary_per_primary" > ${OUT}
    echo -e "${PRIMARY}\t${SECONDARY}\t$(echo "scale=4; ${SECONDARY}/${PRIMARY}" | bc)" >> ${OUT}
}

run_antisense() {
    local BAM=$1 NAME=$2
    local OUT=${METRICS_DIR}/antisense/${NAME}.antisense.txt
    echo "$(date) -- antisense: ${NAME}"
    TOTAL=$(samtools view -@ ${THREADS} -F 0x904 -c ${BAM})
    ANTISENSE=$(samtools view -@ ${THREADS} -F 0x904 -f 0x10 -c ${BAM})
    echo -e "primary_mapped\tantisense\tantisense_rate" > ${OUT}
    echo -e "${TOTAL}\t${ANTISENSE}\t$(echo "scale=4; ${ANTISENSE}/${TOTAL}" | bc)" >> ${OUT}
}

# O(n) awk hash — much faster than sort-based approach for large BAMs
run_aln_per_read() {
    local BAM=$1 NAME=$2
    local OUT=${METRICS_DIR}/aln_per_read/${NAME}.aln_per_read.txt
    echo "$(date) -- aln_per_read: ${NAME}"
    echo -e "n_alignments\tread_count" > ${OUT}
    samtools view -@ ${THREADS} -F 0x804 ${BAM} \
        | awk '
            { count[$1]++ }
            END {
                for (read in count) freq[count[read]]++
                for (n in freq) print n"\t"freq[n]
            }' \
        | sort -k1,1n \
        >> ${OUT}
}

# ════════════════════════════════════════════════════════════════════════
# SECTION 1 — True SR mapping rate from FASTQ read counts
#
# STAR does not write unmapped reads to BAM by default, so flagstat
# shows ~100%. True rate = primary_mapped_BAM / reads_in_FASTQ.
# ════════════════════════════════════════════════════════════════════════
if [ "${RUN_MAPPING_RATE}" -eq 1 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "SECTION 1: True SR mapping rate"
    echo "════════════════════════════════════════════════════════"

    COUNT1=$(awk 'NR%4==1' ${SR_FASTQ_1} | wc -l)
    COUNT2=$(awk 'NR%4==1' ${SR_FASTQ_2} | wc -l)
    TOTAL_INPUT=$((COUNT1 + COUNT2))
    echo "  Sample 1: ${COUNT1} reads"
    echo "  Sample 2: ${COUNT2} reads"
    echo "  Total input: ${TOTAL_INPUT} reads"

    OUT=${METRICS_DIR}/flagstat/SR_fastq_read_counts.txt
    echo -e "reference\ttotal_input_reads\tprimary_mapped\tmapping_rate_pct" > ${OUT}

    for NAME in SR_merged_genome SR_merged_orf_only SR_merged_orf_1000 \
                SR_merged_our_ref SR_merged_nagalakshmi; do
        BAM=${SR_DIR}/${NAME}.bam
        [ ! -f ${BAM} ] && { echo "  WARNING: ${BAM} not found, skipping."; continue; }
        PRIMARY=$(samtools view -@ ${THREADS} -F 0x904 -c ${BAM})
        RATE=$(echo "scale=4; ${PRIMARY} / ${TOTAL_INPUT} * 100" | bc)
        echo -e "${NAME}\t${TOTAL_INPUT}\t${PRIMARY}\t${RATE}" >> ${OUT}
        echo "  ${NAME}: ${PRIMARY} / ${TOTAL_INPUT} = ${RATE}%"
    done
    echo "$(date) -- Section 1 done."
fi

# ════════════════════════════════════════════════════════════════════════
# SECTION 2 — Core alignment metrics for all BAMs
#
# DRS BAMs (minimap2):  multimapping via secondary flag
# SR genome BAM (bwa):  multimapping via NH=0 (always 0; see Section 3)
# SR transcriptome BAMs (STAR): multimapping via NH tag
#
# Note: STAR suppresses antisense alignments by default for transcriptome
# BAMs, so antisense_rate will be 0 or near-0 for those — this is a
# STAR behaviour, not a biological result. Mark as N/A in downstream
# figures (handled by plot_figure3.py / plot_reference_comparison.py).
# ════════════════════════════════════════════════════════════════════════
if [ "${RUN_CORE_METRICS}" -eq 1 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "SECTION 2: Core alignment metrics"
    echo "════════════════════════════════════════════════════════"

    # DRS BAMs — minimap2, secondary-flag multimapping
    DRS_BAMS=(DRS_genome DRS_orf_only DRS_orf_1000 DRS_our_ref DRS_nagalakshmi)
    for NAME in "${DRS_BAMS[@]}"; do
        BAM=${DRS_DIR}/${NAME}.bam
        [ ! -f ${BAM} ] && { echo "WARNING: ${BAM} not found, skipping."; continue; }
        echo "--- ${NAME} ---"
        run_flagstat           ${BAM} ${NAME}
        run_mosdepth           ${BAM} ${NAME}
        run_multimapping_minimap2 ${BAM} ${NAME}
        run_antisense          ${BAM} ${NAME}
        run_aln_per_read       ${BAM} ${NAME}
    done

    # SR genome BAM — bwa-mem2, NH-tag multimapping will give 0
    # (true proxy computed in Section 3 via MAPQ=0)
    echo "--- SR_merged_genome ---"
    BAM=${SR_DIR}/SR_merged_genome.bam
    if [ -f ${BAM} ]; then
        run_flagstat           ${BAM} SR_merged_genome
        run_mosdepth           ${BAM} SR_merged_genome
        run_multimapping_STAR  ${BAM} SR_merged_genome   # will be ~0; overridden by Section 3
        run_antisense          ${BAM} SR_merged_genome
        run_aln_per_read       ${BAM} SR_merged_genome
    fi

    # SR transcriptome BAMs — STAR, NH-tag multimapping
    SR_TXOME_BAMS=(SR_merged_orf_only SR_merged_orf_1000 SR_merged_our_ref SR_merged_nagalakshmi)
    for NAME in "${SR_TXOME_BAMS[@]}"; do
        BAM=${SR_DIR}/${NAME}.bam
        [ ! -f ${BAM} ] && { echo "WARNING: ${BAM} not found, skipping."; continue; }
        echo "--- ${NAME} ---"
        run_flagstat           ${BAM} ${NAME}
        run_mosdepth           ${BAM} ${NAME}
        run_multimapping_STAR  ${BAM} ${NAME}
        run_antisense          ${BAM} ${NAME}
        run_aln_per_read       ${BAM} ${NAME}
    done
    echo "$(date) -- Section 2 done."
fi

# ════════════════════════════════════════════════════════════════════════
# SECTION 3 — MAPQ=0 multi-mapping proxy for genome BAMs
#
# bwa-mem2 does NOT write secondary alignment records (flag 0x100)
# regardless of any flags, so the NH/secondary-based metrics in
# Section 2 always give 0 for genome BAMs. The correct proxy is:
#   multi-mapping reads ≈ primary reads with MAPQ = 0
# (bwa-mem2 assigns MAPQ=0 when a read aligns equally well to ≥2 loci)
#
# This section writes <stem>.bwa_multimapping.txt for the genome BAMs.
# These files are read preferentially by plot_reference_comparison.py
# / plot_figure3.py and override the NH-based rate for Genome entries.
# ════════════════════════════════════════════════════════════════════════
if [ "${RUN_BWA_MAPQ0}" -eq 1 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "SECTION 3: MAPQ=0 multi-mapping proxy for genome BAMs"
    echo "════════════════════════════════════════════════════════"

    mapq0_proxy() {
        local BAM=$1 NAME=$2 METHOD=$3
        local OUT=${METRICS_DIR}/flagstat/${NAME}.bwa_multimapping.txt
        echo "$(date) -- MAPQ=0 proxy: ${NAME}"
        TOTAL=$(samtools view -@ ${THREADS} -F 0x904 -c ${BAM})
        MAPQ_GE1=$(samtools view -@ ${THREADS} -F 0x904 -q 1 -c ${BAM})
        MULTI=$((TOTAL - MAPQ_GE1))
        RATE=$(echo "scale=4; ${MULTI} / ${TOTAL}" | bc)
        echo "  MAPQ=0: ${MULTI} / ${TOTAL} = ${RATE}"
        echo -e "total_primary\tmapq0_reads\tmapq0_rate\tmethod" > ${OUT}
        echo -e "${TOTAL}\t${MULTI}\t${RATE}\t${METHOD}" >> ${OUT}
    }

    BAM=${SR_DIR}/SR_merged_genome.bam
    [ -f ${BAM} ] && mapq0_proxy ${BAM} SR_merged_genome bwa_mapq0_proxy

    BAM=${DRS_DIR}/DRS_genome.bam
    [ -f ${BAM} ] && mapq0_proxy ${BAM} DRS_genome minimap2_mapq0_proxy

    echo "$(date) -- Section 3 done."
fi

# ════════════════════════════════════════════════════════════════════════
# SECTION 4 — Soft-clip statistics (DRS only)
#
# Calls softclip_metrics.py. Not meaningful for short reads.
# Output: one TSV per DRS BAM in ${METRICS_DIR}/softclip/
# ════════════════════════════════════════════════════════════════════════
if [ "${RUN_SOFTCLIP}" -eq 1 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "SECTION 4: Soft-clip metrics (DRS only)"
    echo "════════════════════════════════════════════════════════"

    for NAME in DRS_genome DRS_orf_only DRS_orf_1000 DRS_our_ref DRS_nagalakshmi; do
        BAM=${DRS_DIR}/${NAME}.bam
        [ ! -f ${BAM} ] && { echo "WARNING: ${BAM} not found, skipping."; continue; }
        echo "$(date) -- softclip_metrics: ${NAME}"
        python3 ${SCRIPTS_DIR}/softclip_metrics.py \
            --bam       ${BAM} \
            --out       ${METRICS_DIR}/softclip/${NAME}.softclip.tsv \
            --min-reads 5
    done
    echo "$(date) -- Section 4 done."
fi

# ════════════════════════════════════════════════════════════════════════
# SECTION 5 — Metagene coverage profiles
#
# Calls compute_metagene.py for all 4 references × 2 data types × 2
# anchor points (TSS / TES) = 16 profile TSVs.
#
# All references are restricted to the same shared gene set (genes
# with valid UTR data in both our annotation AND the reference UTR
# table) for a fair cross-reference comparison.
#
# Contig naming conventions:
#   orf_only / orf_1000 : bare gene names  (--name_suffix "")
#   our_ref / nagalakshmi : gene_mRNA names (--name_suffix _mRNA)
# ════════════════════════════════════════════════════════════════════════
if [ "${RUN_METAGENE}" -eq 1 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "SECTION 5: Metagene coverage profiles"
    echo "════════════════════════════════════════════════════════"

    # Compute shared gene list (genes with valid UTR data in both sources)
    SHARED_GENES=${PROFILES_DIR}/shared_genes.txt
    echo "$(date) -- Computing shared gene list..."
    python3 - << PYEOF
import pandas as pd

our = pd.read_csv("${OUR_UTR}", sep="\t")
our.columns = ["Gene"] + list(our.columns[1:])

nag = pd.read_csv("${NAGALAKSHMI_UTR}")
nag.columns = ["Gene"] + list(nag.columns[1:])
nag_utr5 = pd.to_numeric(nag.iloc[:, 1], errors="coerce")
nag_utr3 = pd.to_numeric(nag.iloc[:, 2], errors="coerce")
nag["_u5"] = nag_utr5
nag["_u3"] = nag_utr3

nag_valid  = set(nag.dropna(subset=["_u5","_u3"])["Gene"])
our_genes  = set(our["Gene"])
shared     = sorted(nag_valid & our_genes)

with open("${SHARED_GENES}", "w") as fh:
    fh.write("\n".join(shared) + "\n")

print(f"Shared gene set: {len(shared)} genes")
print(f"  Our annotation total : {len(our_genes)}")
print(f"  Nagalakshmi valid UTR: {len(nag_valid)}")
print(f"  Only in ours (unused): {len(our_genes - nag_valid)}")
PYEOF
    echo "$(date) -- Shared genes written: ${SHARED_GENES}"

    # Helper: run compute_metagene.py for one BAM
    run_metagene() {
        local BAM=$1
        shift
        echo ""
        echo "$(date) -- metagene: $(basename ${BAM} .bam)"
        python3 ${SCRIPTS_DIR}/compute_metagene.py \
            --bam         ${BAM} \
            --outdir      ${PROFILES_DIR} \
            --window      1000 \
            --min_orf     200 \
            --restrict_to ${SHARED_GENES} \
            "$@"
    }

    # Short Read
    echo ""; echo "=== Short Read ==="
    run_metagene ${SR_DIR}/SR_merged_orf_only.bam \
        --name_suffix ""

    run_metagene ${SR_DIR}/SR_merged_orf_1000.bam \
        --fixed_utr5 1000 --fixed_utr3 1000 \
        --name_suffix ""

    run_metagene ${SR_DIR}/SR_merged_our_ref.bam \
        --utr_file ${OUR_UTR} \
        --utr5_col final_five_prime_utr --utr3_col final_three_prime_utr \
        --name_suffix _mRNA

    run_metagene ${SR_DIR}/SR_merged_nagalakshmi.bam \
        --utr_file ${NAGALAKSHMI_UTR} \
        --utr5_col five_prime_utr --utr3_col three_prime_utr \
        --name_suffix _mRNA

    # DRS
    echo ""; echo "=== DRS ==="
    run_metagene ${DRS_DIR}/DRS_orf_only.bam \
        --name_suffix ""

    run_metagene ${DRS_DIR}/DRS_orf_1000.bam \
        --fixed_utr5 1000 --fixed_utr3 1000 \
        --name_suffix ""

    run_metagene ${DRS_DIR}/DRS_our_ref.bam \
        --utr_file ${OUR_UTR} \
        --utr5_col final_five_prime_utr --utr3_col final_three_prime_utr \
        --name_suffix _mRNA

    run_metagene ${DRS_DIR}/DRS_nagalakshmi.bam \
        --utr_file ${NAGALAKSHMI_UTR} \
        --utr5_col five_prime_utr --utr3_col three_prime_utr \
        --name_suffix _mRNA

    echo "$(date) -- Section 5 done."
    echo "$(date) -- Download ${PROFILES_DIR}/*.tsv and run plot_figure4.py locally."
fi

echo ""
echo "════════════════════════════════════════════════════════"
echo "$(date) -- All sections complete."
echo "Metrics directory: ${METRICS_DIR}"
echo "════════════════════════════════════════════════════════"
