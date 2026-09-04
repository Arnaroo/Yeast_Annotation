#!/usr/bin/env bash
#
# read_migration.sh
# ==================
# Follows individual reads across the two references, to test directly
# whether an extended UTR reaching into a neighbour's territory takes
# reads that belong to the neighbour. A gene-level abundance correlation
# (abundance_concordance.py) answers this only indirectly: two genes
# swapping reads symmetrically would still sit on the diagonal.
#
# Every read is aligned to both references in the same run with the same
# settings, and its best-hit transcript is recorded under each. A read
# assigned to gene A under the prior reference and to gene B under the
# released one has migrated; read_migration_analyse.py classifies every
# migration and splits out the subset attributable to the newly
# overlapping gene pairs found by neighbour_overlap.py.
#
# Memory is kept bounded by sorting on disk rather than joining in
# memory -- the short-read set can be tens of millions of reads.
#
# Usage:
#   read_migration.sh <ref_new.fa> <ref_old.fa> <drs.fq.gz> <sr.fq.gz> <outdir>
#
# Output: <outdir>/drs.transitions.tsv, <outdir>/sr.transitions.tsv
#         two columns joined into: old_transcript  new_transcript  n_reads

set -euo pipefail

REF_NEW="$1"
REF_OLD="$2"
DRS_FQ="$3"
SR_FQ="$4"
OUT="$5"

THREADS="${THREADS:-10}"
SORTMEM="${SORTMEM:-2G}"
SORTTMP="$OUT/sorttmp"
mkdir -p "$OUT" "$SORTTMP"

echo "read_migration.sh"
echo "started $(date -Is)"
echo

# assign <label> <ref> <fq> <preset>
#   emits  readname <tab> transcript_or_UNMAPPED, sorted by read name
assign () {
    local label="$1" ref="$2" fq="$3" preset="$4"
    local out="$OUT/${label}.assign.tsv.gz"
    if [ -s "$out" ]; then
        echo "skip $label, assignments present"
        return
    fi
    echo "assigning $label"
    local t0=$SECONDS
    minimap2 -ax "$preset" -k14 -N10 --secondary=yes -t "$THREADS" \
             "$ref" "$fq" 2> "$OUT/${label}.assign.mm2.err" \
      | samtools view -@ 2 -F 0x900 - \
      | awk -v OFS='\t' '{
            unmapped = int($2 / 4) % 2
            print $1, (unmapped ? "UNMAPPED" : $3)
        }' \
      | sort -S "$SORTMEM" -T "$SORTTMP" -k1,1 \
      | gzip -1 > "$out"
    echo "  done in $((SECONDS - t0))s, $(zcat "$out" | wc -l) reads"
}

assign drs_new "$REF_NEW" "$DRS_FQ" map-ont
assign drs_old "$REF_OLD" "$DRS_FQ" map-ont
assign sr_new  "$REF_NEW" "$SR_FQ"  sr
assign sr_old  "$REF_OLD" "$SR_FQ"  sr

# ------------------------------------------------------------------ join
# Paste the two assignments per technology into one table of
#   old_transcript  new_transcript  n_reads
# small enough to analyse in read_migration_analyse.py.
for tech in drs sr; do
    out="$OUT/${tech}.transitions.tsv"
    if [ -s "$out" ]; then
        echo "skip $tech transitions, present"
        continue
    fi
    echo "joining $tech"
    join -t $'\t' -j 1 \
        <(zcat "$OUT/${tech}_old.assign.tsv.gz") \
        <(zcat "$OUT/${tech}_new.assign.tsv.gz") \
      | awk -v OFS='\t' '{ print $2, $3 }' \
      | sort -S "$SORTMEM" -T "$SORTTMP" \
      | uniq -c \
      | awk -v OFS='\t' 'BEGIN { print "old","new","n_reads" } { print $2, $3, $1 }' \
      > "$out"
    echo "  $(( $(wc -l < "$out") - 1 )) distinct transitions"
done

rm -rf "$SORTTMP"
echo
echo "finished $(date -Is)"
