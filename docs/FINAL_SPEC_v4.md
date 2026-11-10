# BLAST validation v4: finalized spec for the rerun (2026-11-09)

Decided with Amanda on 2026-11-09. Supersedes the name-based read-level method (v3) for the rerun.

## Why a rerun
The v2/v3 BLAST TSVs have `sscinames = N/A` on every hit: BLAST could not find `taxdb` (it sits next to nt,
but `$BLASTDB` was never set). Hit names therefore came from free-text titles, and 44% of MISASSIGNED reads
in the top 30 taxa were the right species under another name (Klebsiella ornithinolytica for Raoultella,
"[Mycobacterium] chelonae", "MAG: ..."). NCBI itself now files taxids 54291 / 575 under Klebsiella, so name
matching against Kraken's older names cannot be fixed with names alone. v2.1 consequence: both Raoultella
species had 0 TRUE_POSITIVE detections and are dropped from the Time-Series primary filter.

## Decisions
| Item | Decision |
|---|---|
| Matching | **NCBI taxid rolled up to species** (`staxids` in outfmt; taxonomy = taxdmp_2026-09-01, matching nt 2026_09). Genus-level Kraken targets match at the genus. Merged taxids resolved. |
| Same-complex sibling wins | **EXCLUDED_COMPLEX** (not MISASSIGNED). Complex = NCBI "species group"/"species subgroup" (broadest below genus). Only when *every* top-bitscore named hit is a sibling. |
| Near-tie margin | **Decide after the rerun**: sweep 2% / 5% / 10% on the corrected output; Hannah picks. |
| Paper input | **Read-level supported detections + corrected relative abundance** replace the v2.1 TRUE_POSITIVE rule in Time_Series `config/ww_filter.yml`, after Amanda reviews old-vs-new differences. |
| Kept as is | nt 2026_09; up to 500 reads per read set (seqtk seed 42, same reads as v2); max_target_seqs 100; %id >= 80; supported = CONFIRMED > 2; uncultured/environmental hits ignored. |

## Synonym table (added 2026-11-09, Amanda: "we really need a synonym table")
`resources/taxon_synonyms.csv` (built by `tools/build_synonym_table.py` from taxdmp_2026-09-01 + `resources/taxon_synonyms_manual.csv`):
one row per (target taxid, name): Kraken name, NCBI current name, NCBI synonyms / equivalent names, and the
species epithet under each NCBI synonym of its genus. 964 names for 142 targets; NCBI has renamed 2 targets
(Raoultella ornithinolytica / planticola -> Klebsiella). Safety rules (each found on real data):
- Subspecies-form synonyms are dropped ("Aeromonas hydrophila subsp. anaerogenes" for A. caviae would make
  A. hydrophila match A. caviae).
- A synonym that is now the valid name of ANOTHER cellular species is dropped (NCBI keeps "Klebsiella aerogenes"
  under K. pneumoniae, "Mycobacterium tuberculosis" under M. avium, "Proteus vulgaris" under P. penneri).
  Vectors and host-named viruses are ignored for this check. Dropped names: `resources/taxon_synonyms_dropped.csv`.
- Title prefixes "MAG:", "TPA_asm:", "TPA_exp:" are stripped before matching.
Used by every name match (`--synonyms`): the hit-without-taxid fallback in taxid mode, and name mode on the
existing v3 TSVs. On the v3 misassigned reads (top 30 taxa) it recovers 40,690 of 92,300 (44.1%): all 40,539
name-only cases plus Peptoclostridium/Clostridium difficile, Enterobacter aerogenes, Aeromonas punctata.
Edit `taxon_synonyms_manual.csv` (taxid,synonym,note) to add names; rebuild with the tool.

## Open (needs Amanda/Hannah)
- **Denominator for the validated fraction.** v3 rule: CONFIRMED / all reads with a hit (EXCLUDED counts against
  the taxon). Both are reported: `validated_fraction` (v3 rule) and `validated_fraction_informative` =
  CONFIRMED / (CONFIRMED + MISASSIGNED), which treats EXCLUDED and EXCLUDED_COMPLEX as uninformative. Labelling a
  sibling read EXCLUDED only changes abundance under the second definition.

## Per-read rule (v4)
1. Drop uncultured/environmental hits (by taxonomy name / "environmental samples" lineage).
2. No named hit left -> **EXCLUDED**.
3. Target-species hit with bitscore >= (1 - margin) x top named bitscore and %id >= 80 -> **CONFIRMED**.
4. Otherwise, all top-bitscore named hits are sibling species in the target's NCBI complex -> **EXCLUDED_COMPLEX**.
5. Otherwise -> **MISASSIGNED**.
Per detection: corrected_abundance = validated_fraction x Kraken clade reads; relative = / total QC'd dehosted reads.

## Code (implemented 2026-11-09)
- `workflow/scripts/blast_sample.sh`: `export BLASTDB=<nt dir>`, outfmt adds `staxids`; fails fast if taxdb missing.
- `workflow/scripts/taxonomy.py`: species / complex / merged / uncultured lookups.
- `workflow/scripts/false_positive_detection.py`: `--taxonomy` turns on taxid mode (`per_read_frame_tax`,
  `label_reads`, `summarize_labels`); reads 7-col (v3) or 8-col (v4) TSVs; new columns `rl_excluded_complex`,
  `validated_fraction_informative`. Legacy `classify()` untouched.
- `tools/near_tie_margin_sweep.py --taxonomy`: same per-read code, all margins in one pass.
- `tools/build_taxonomy_table.py`: taxdump -> `/n/netscratch/hhealy_lab/avdarling/ncbi_taxonomy_2026-09-01/taxonomy_table.tsv`.
- `tests/test_taxid_mode.py`: 8 tests (rename, MAG, merged, strain, sibling, mixed tie, uncultured, low id).
- `tests/test_synonyms.py`: 4 tests (renames, prefixes, no leaks across species incl. A. hydrophila/A. caviae).

## Rerun plan
1. Pilot: one sample (GKWWBAA00014) -> check sscinames + staxids populated, runtime.
2. Full BLAST rerun from the existing extracted FASTAs to a NEW results dir on netscratch (v2 left intact).
3. Classify + 2/5/10% sweep with `--taxonomy` and `--clade_counts`.
4. Old-vs-new comparison per taxon (supported detections, abundance), for Amanda before anything reaches the paper.
