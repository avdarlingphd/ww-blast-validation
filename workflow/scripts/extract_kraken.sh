#!/bin/bash
# Extract the reads Kraken2 assigned to one (sample, taxid) into per-read-set FASTAs (kraken variant).
#
# Uses KrakenTools extract_kraken_reads.py with --include-children so reads classified to any
# strain/subspecies under TAXID are captured (the .kreport is walked as a taxonomy tree). The
# extracted FASTQs are converted to FASTA with seqtk.
#
# Zero-read handling: if no reads are classified to TAXID (or KrakenTools fails / the taxid is
# absent from the report), BOTH output FASTAs are still created EMPTY. Downstream, blast_pair emits
# an empty result and the classifier reads it as NO_DATA rather than crashing. This is exactly the
# case that produced the original 1,415 NO_DATA calls.
#
# Env:  PYTHON  KRAKENTOOLS  SEQTK
# Args: SAMPLE TAXID KRAKEN KREPORT FASTQ_R1 FASTQ_R2 OUT_R1_FASTA OUT_R2_FASTA
set -uo pipefail   # NOT -e: a failed/empty extraction must still leave empty outputs, not abort
SAMPLE="$1"; TAXID="$2"; KRAKEN="$3"; KREPORT="$4"; FQ1="$5"; FQ2="$6"; OUT_R1="$7"; OUT_R2="$8"
: "${PYTHON:?}"; : "${KRAKENTOOLS:?}"; : "${SEQTK:?}"

od="$(dirname "$OUT_R1")"; mkdir -p "$od"
tmp1="${od}/${SAMPLE}_taxid_${TAXID}_R1.fastq"
tmp2="${od}/${SAMPLE}_taxid_${TAXID}_R2.fastq"

echo "extract: sample=${SAMPLE} taxid=${TAXID} (include-children)"
"$PYTHON" "$KRAKENTOOLS" \
    -k "$KRAKEN" -r "$KREPORT" \
    -s "$FQ1" -s2 "$FQ2" \
    -t "$TAXID" --include-children --fastq-output \
    -o "$tmp1" -o2 "$tmp2" || echo "  KrakenTools returned non-zero for taxid ${TAXID}; treating as 0 reads"

# FASTQ -> FASTA (empty in, empty out). Always leave both outputs present.
if [ -s "$tmp1" ]; then "$SEQTK" seq -a "$tmp1" > "$OUT_R1"; else : > "$OUT_R1"; fi
if [ -s "$tmp2" ]; then "$SEQTK" seq -a "$tmp2" > "$OUT_R2"; else : > "$OUT_R2"; fi
rm -f "$tmp1" "$tmp2"

n1=$(grep -c "^>" "$OUT_R1" 2>/dev/null || echo 0)
echo "extract done: ${SAMPLE} taxid ${TAXID} -> ${n1} reads"
