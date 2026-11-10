# Changelog

Version history for the pipeline. The README describes current behavior; this file records how it
got there.

## v4 (2026-11-09): taxid matching, synonym table, sibling-complex label
- **BLAST now sees NCBI taxonomy.** `blast_sample.sh` exports `BLASTDB` to the nt directory (where `taxdb.btd/.bti`
  live) and fails fast if taxdb is missing. Before this, `sscinames` was `N/A` for every hit, so hit names came
  from free-text titles ("Klebsiella ornithinolytica", "[Mycobacterium] chelonae", "MAG: ..."), and correct
  reads were scored MISASSIGNED. The outfmt also adds `staxids` (now 8 columns; 7-column v3 files still load).
- **Taxid matching** (`--taxonomy`, `workflow/scripts/taxonomy.py`): a hit matches the target if its NCBI taxid
  rolls up to the target species (merged taxids resolved; genus-level targets match at the genus). Taxonomy
  table built by `tools/build_taxonomy_table.py` from a dated taxdump matching the nt snapshot.
- **Synonym table** (`--synonyms`, `resources/taxon_synonyms.csv`, `tools/build_synonym_table.py`): every name a
  target goes by, from NCBI plus `taxon_synonyms_manual.csv`. Drops subspecies-form synonyms and synonyms that
  are now another species' valid name (`taxon_synonyms_dropped.csv`). "MAG:" / "TPA_asm:" title prefixes are
  stripped before any name match.
- **EXCLUDED_COMPLEX** read label: not confirmed, and every top-scoring named hit is a sibling species in the
  target's NCBI complex ("species group"). New report columns `rl_excluded_complex`,
  `validated_fraction_informative` (CONFIRMED / (CONFIRMED + MISASSIGNED)).
- `tools/near_tie_margin_sweep.py` shares the per-read code (`per_read_frame_tax`, `label_reads`).
- Tests: `tests/test_taxid_mode.py` (8), `tests/test_synonyms.py` (4). Decisions: `docs/DECISION_POINTS.md`;
  spec: `docs/FINAL_SPEC_v4.md`; figure: `docs/master_methods_figure.png`.

## v3: read-level validation (additive)
- Each read is CONFIRMED (target within a near-tie bitscore margin of the top named hit, %id >= 80;
  uncultured ignored), MISASSIGNED, or EXCLUDED (no named hit). Per detection: `validated_fraction`,
  `supported` (CONFIRMED > 2), and with `--clade_counts` a corrected / relative abundance
  (`tools/build_clade_counts.py`). `tools/near_tie_margin_sweep.py` runs several margins in one pass.
  The legacy per-detection classification is unchanged.

## v2.3
- Optional **tie-aware matching** (`tie_aware` / `--tie_aware`): a read counts as a match if the target
  is among its tied-best (max-bitscore) hits, instead of only the single arbitrary top hit. Removes the
  database-over-representation bias that lets a more-deposited sibling win a tie (e.g. *K. pneumoniae* over
  a true *K. variicola*). Reads that match a different organism strictly better are still excluded.
- Optional **uncultured-as-uninformative** (`uncultured_uninformative` / flag): reads whose best hit is an
  uncultured/environmental entry (and not a target match) are dropped from the denominator rather than
  counted as non-matches. UNCULTURED_DOMINANT remains as a separate flag.
- Both default **off** (backward-compatible); recorded in PROVENANCE.txt (`matching_modes`). 4 new unit tests.

## v2.2
- Optional **species-complex scoring**: a user-supplied `complex_crosswalk` (complex_name, member_organism)
  adds complex-level columns (`complex_name`, `complex_classification`, `avg_complex_pct_match`) alongside
  the species-level call, so complex members split across siblings (ACB, E. cloacae, K. pneumoniae, …) are
  recognized at the complex level. Name-based, backward-compatible (blank crosswalk = species-level only).
  Ships a PROVISIONAL `resources/complex_crosswalk.example.csv` to be replaced by the real
  RPIP-panel / clinical complex scheme.
- Added `config/config.example.yaml` template (placeholder paths) to make the tool easy to reuse elsewhere.

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
