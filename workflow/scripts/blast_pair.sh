#!/bin/bash
# BLAST one (sample, taxid) read pair. Subsample to MAX_READS, then blastn vs nt with -max_target_seqs
# MAX_TARGET and bitscore in the outfmt (needed for top-hit-by-bitscore classification).
# Env: BLASTN SEQTK DB [MAX_READS=500] [MAX_TARGET=100] [THREADS=8]
# Args: IN_R1 IN_R2 OUT_R1 OUT_R2
set -euo pipefail
IN_R1="$1"; IN_R2="$2"; OUT_R1="$3"; OUT_R2="$4"
: "${BLASTN:?}"; : "${SEQTK:?}"; : "${DB:?}"
MAX_READS="${MAX_READS:-500}"; MAX_TARGET="${MAX_TARGET:-100}"; THREADS="${THREADS:-8}"
OUTFMT="6 qseqid sscinames pident length evalue bitscore stitle"
od="$(dirname "$OUT_R1")"; mkdir -p "$od"

run_one () {
    local infa="$1" out="$2" q="$1" n
    # Empty/near-empty FASTA -> emit an empty result (classifier reads it as NO_DATA rather than failing)
    n=$(grep -c "^>" "$infa" 2>/dev/null || echo 0)
    if [ "${n:-0}" -eq 0 ]; then
        : > "$out"; echo "  ${infa}: 0 reads -> empty result"; return 0
    fi
    if [ "${n}" -gt "${MAX_READS}" ]; then
        q="${od}/$(basename "$out" .tsv)_sub.fasta"
        "$SEQTK" sample -s 42 "$infa" "$MAX_READS" > "$q"
    fi
    # -task megablast (word size 28): the fast, high-identity mode, appropriate for validating short
    # Illumina reads against nt (10-50x faster than -task blastn). Use dc-megablast/blastn only if you
    # need to catch divergent hits, which read-level taxonomic validation does not.
    "$BLASTN" -query "$q" -db "$DB" -task megablast -num_threads "$THREADS" \
        -max_target_seqs "$MAX_TARGET" -outfmt "$OUTFMT" -out "$out"
}
run_one "$IN_R1" "$OUT_R1"
run_one "$IN_R2" "$OUT_R2"
echo "done: $OUT_R1 / $OUT_R2"
