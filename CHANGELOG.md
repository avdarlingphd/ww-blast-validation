# Changelog

Version history for the pipeline. The README describes current behavior; this file records how it
got there.

## v2.1
- **Species-level organism matching** in the classifier: a species-level candidate requires the hit to
  agree at BOTH genus and species epithet; a genus-only or sibling-species hit no longer confirms it.
- **`TRUE_POSITIVE` requires > 2 matching reads per set** (`min_match_reads: 3`); the v2 single-read
  "rescue" (`min_match_reads: 1`) was reverted. A single-read call is never positive (read floor).
- **BLAST runs once per sample** (`blast_sample`), tagging each read by its taxid, instead of once per
  `(sample, taxid)` pair — the large `nt` load happens ~once-per-sample instead of ~once-per-pair.
- **`-task megablast`** (was `-task blastn`) for fast, high-identity read validation.
- **Automatic per-run provenance** (`PROVENANCE.txt` + `config.snapshot.yaml`); `nt` pinnable to a dated snapshot.
- **Unit tests** for the classification logic (`tests/`).

## v2
- Top-hit-per-read-by-bitscore classification (bitscore added to the BLAST outfmt).
- `-max_target_seqs 100` (was 20) so the true organism is not crowded out of the hit list.
- Read extraction folded into the Snakemake workflow (previously pre-extracted FASTAs on scratch).
- Per-sample batched extraction: read the large `.kraken` once per sample, not once per taxid.

## v1
- Original per-`(sample, taxid)` shell scripts (no workflow manager).

## Internal note (Healy lab)
- Older loose copies of `false_positive_detection.py` under `.../Linux Scripts/` and
  `.../Biofilm-Project/4_blastn_validation/` are deprecated v1 and must not be run — this repo is canonical.
