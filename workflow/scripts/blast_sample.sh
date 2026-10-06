#!/bin/bash
# BLAST all candidate reads for ONE sample in a single pass, so the ~1 TB nt database is memory-mapped
# ONCE per sample (~250x over the run) instead of once per (sample, taxid) pair (~8,600x). That DB load
# is the dominant cost, so batching per sample is the big speedup.
#
# For each of the sample's candidate taxids, each read set (R1, R2) is subsampled to MAX_READS and its
# read headers are TAGGED with the taxid and direction as "<taxid>|<R1|R2>|<original_read_id>". All
# tagged reads are concatenated into one query and BLASTed once. false_positive_detection.py then splits
# the single per-sample result back out by that tag and classifies each (sample, taxid) as before.
#
# Env:  BLASTN SEQTK DB [MAX_READS=500] [MAX_TARGET=100] [THREADS=8]
# Args: SAMPLE  FASTA_TREE  OUT_TSV
set -euo pipefail
SAMPLE="$1"; FASTA_TREE="$2"; OUT="$3"
: "${BLASTN:?}"; : "${SEQTK:?}"; : "${DB:?}"
MAX_READS="${MAX_READS:-500}"; MAX_TARGET="${MAX_TARGET:-100}"; THREADS="${THREADS:-8}"
OUTFMT="6 qseqid sscinames pident length evalue bitscore stitle"
od="$(dirname "$OUT")"; mkdir -p "$od"

tmpq="$(mktemp "${od}/.${SAMPLE}_query.XXXXXX.fasta")"
trap 'rm -f "$tmpq"' EXIT

# Tag every read header with "<taxid>|<read>|<id>" (first whitespace token of the original id).
tag () {  # args: taxid read   (reads FASTA on stdin)
    awk -v t="$1" -v r="$2" '/^>/{print ">" t "|" r "|" substr($1,2); next} {print}'
}

for d in "${FASTA_TREE}/blast_validation_${SAMPLE}_"*/ ; do
    [ -d "$d" ] || continue
    taxid="${d%/}"; taxid="${taxid##*_}"                 # trailing _<taxid> of the dir name
    for read in R1 R2; do
        fa="${d}${SAMPLE}_taxid_${taxid}_${read}.fasta"
        [ -s "$fa" ] || continue
        n=$(grep -c "^>" "$fa" 2>/dev/null || echo 0)
        [ "${n:-0}" -gt 0 ] || continue
        if [ "$n" -gt "$MAX_READS" ]; then
            "$SEQTK" sample -s 42 "$fa" "$MAX_READS" | tag "$taxid" "$read" >> "$tmpq"
        else
            tag "$taxid" "$read" < "$fa" >> "$tmpq"
        fi
    done
done

# No candidate reads for this sample -> empty result (classify reads it as NO_DATA for every taxid).
if [ ! -s "$tmpq" ]; then
    : > "$OUT"; echo "  ${SAMPLE}: no candidate reads -> empty result"; exit 0
fi

"$BLASTN" -query "$tmpq" -db "$DB" -task megablast -num_threads "$THREADS" \
    -max_target_seqs "$MAX_TARGET" -outfmt "$OUTFMT" -out "$OUT"
echo "done: $OUT ($(grep -c '^>' "$tmpq") tagged query reads)"
