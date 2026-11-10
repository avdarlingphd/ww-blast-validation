# BLAST read-level validation: decision points (v4, 2026-11-09)

One page for the advisor meeting. Each row is a choice in the method, the options considered, what we chose
and why. Methods figure: `docs/master_methods_figure.png`. Full spec: `docs/FINAL_SPEC_v4.md`.

## A. Decided

| # | Decision point | Options considered | Chosen | Why |
|---|---|---|---|---|
| 1 | Unit of validation | per detection (TRUE/FALSE POSITIVE buckets, v2.1) vs per read | **Per read**, then re-aggregate | A detection is a mix of right and wrong reads; per-read gives a fraction that corrects abundance instead of an all-or-nothing call. |
| 2 | Reads checked per detection | all reads vs a subsample | **Up to 500 per read set** (seqtk seed 42) | Subsample only estimates the confirmed fraction; no reads are discarded from abundance. Same design as NAO/SecureBio mgs-workflow. |
| 3 | Reference database | nt vs curated genomes | **NCBI nt, dated snapshot 2026_09** | Broadest reference, so a wrong Kraken call has the chance to lose to the right organism. Pinned for reproducibility. |
| 4 | Hits kept per read | 20 vs 100 | **100** (`max_target_seqs`) | The cap is a collection heuristic (BLAST keeps the first N found, then sorts), so a small cap can crowd the true organism out of the list. |
| 5 | Ranking hits | raw % identity vs bitscore | **Bitscore** | Bitscore combines identity and alignment length, so a long 98% alignment beats a short 100% one; within a read it orders hits the same way as E-value. |
| 6 | Uncultured / environmental hits | count against the target vs ignore | **Ignore** | They reflect database composition, not a competing named organism. |
| 7 | How a hit is matched to the target | names from titles vs NCBI taxonomy names vs taxids | **NCBI taxid rolled up to species**, with a synonym table as fallback | Names are unstable: NCBI now files Raoultella under Klebsiella with the same taxids. Taxids cannot drift; the synonym table covers hits without a taxid. |
| 8 | Synonym table rules | all NCBI synonyms vs filtered | **Filtered**: drop subspecies-form synonyms and any synonym that is now another species' valid name | Unfiltered, A. hydrophila would match A. caviae and K. aerogenes would match K. pneumoniae. 24 names dropped, listed in `resources/taxon_synonyms_dropped.csv`. |
| 9 | Minimum identity | none / 80% / higher | **80%** | Below this a "match" is not species-level evidence. |
| 10 | Sibling in the same complex wins | MISASSIGNED vs EXCLUDED vs count as target | **EXCLUDED_COMPLEX** (complex = NCBI "species group") | Species inside a complex are often indistinguishable on 150 bp; a sibling win is not evidence against the taxon. Uses NCBI's own complex definitions, not a hand list. |
| 11 | Corrected abundance | confirmed reads only vs fraction x Kraken count | **validated fraction x Kraken clade reads / total QC'd dehosted reads** | Keeps Kraken's full depth; BLAST only supplies the correction factor. |
| 12 | Supported detection | >= 1 / > 2 confirmed reads | **> 2 confirmed reads** | One or two reads are too easily chance hits. |
| 13 | What feeds the Time-Series paper | v2.1 TRUE_POSITIVE rule vs read-level | **Read-level supported detections + corrected abundance** | Applied only after Amanda reviews old vs new per taxon. |

## B. Still open

| # | Decision point | Options | Status |
|---|---|---|---|
| 14 | Near-tie margin | 2% / 5% / 10% (NAO uses 10%) | Rerun computes all three; Hannah picks. On the uncorrected data, 2% vs 10% changed supported detections by 0.8%. |
| 15 | Denominator of the validated fraction | CONFIRMED / all reads with a hit vs CONFIRMED / (CONFIRMED + MISASSIGNED) | Both are computed. Only the second lets EXCLUDED / EXCLUDED_COMPLEX reads stop counting against the taxon. |

## C. Quality checks and what they caught

Each check below found something; each is fixed and covered by tests (35 tests pass).

| Check | Finding | Fix |
|---|---|---|
| Which organism wins the MISASSIGNED reads (`tools/misassigned_winners.py`) | BLAST could not see NCBI's taxonomy file, so every hit was named from its free-text title. 44% of misassigned reads (top 30 taxa) were the correct species under another name (Raoultella/Klebsiella, "[Mycobacterium] chelonae", "MAG:" titles). | BLAST now finds taxdb and also reports taxids (#7); synonym table (#8). |
| Synonym table audit against all NCBI species | Old synonyms that are now other species' names would merge distinct species. | Filter rules (#8). |
| Pilot rerun, sample GKWWBAA00014 (39 detections, 19,605 reads) | Every hit now has an NCBI name and taxid. Old BLAST + synonym table and new BLAST + taxids give **identical CONFIRMED counts** at 2/5/10% (19,348 / 19,418 / 19,463). MISASSIGNED at 2%: 1,832 before the fix, 226 after (+15 EXCLUDED_COMPLEX). | Two independent routes agree, so the fix is confirmed. |

## D. Next steps
1. Full BLAST rerun, 250 samples (~4.6 h per sample, run in parallel) to a new netscratch folder; v2 kept intact.
2. Classify + 2/5/10% sweep on the new output.
3. Old vs new per taxon (supported detections, abundance) for Amanda, then Hannah picks #14 and #15.
4. Then update the Time-Series wastewater filter (#13).
