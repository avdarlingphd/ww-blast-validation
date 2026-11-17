# BLAST read-level validation: decision points (v4, 2026-11-09; updated 2026-11-17)

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
| 9 | Minimum identity | none / 80% / 90% | **90%** (raised from 80%, Amanda 2026-11-17) | A confirmed C. difficile read rested on an 87.8% match to a plasmid; below 90% a match is not species-level evidence. Most confirmed reads are >= 97%. |
| 10 | Sibling in the same complex wins | MISASSIGNED vs EXCLUDED vs count as target | **EXCLUDED_COMPLEX** (complex = NCBI "species group") | Species inside a complex are often indistinguishable on 150 bp; a sibling win is not evidence against the taxon. Uses NCBI's own complex definitions, not a hand list. |
| 11 | Corrected abundance | confirmed reads only vs fraction x Kraken count | **validated fraction x Kraken clade reads / total QC'd dehosted reads** | Keeps Kraken's full depth; BLAST only supplies the correction factor. |
| 12 | Minimum confirmed reads per detection | in the validation vs downstream; 1 vs 3 | **Not part of the validation.** The report gives CONFIRMED counts; a minimum (default 3) is an option in the Time-Series R script that can be switched off (Amanda 2026-11-17). A relative-abundance threshold instead was considered and not adopted. | Keeps validation strictly read by read. |
| 13 | What feeds the Time-Series paper | v2.1 TRUE_POSITIVE rule vs read-level | **Read-level supported detections + corrected abundance** | Applied only after Amanda reviews old vs new per taxon. |
| 14 | Near-tie margin | 2% / 5% / 10% (NAO uses 10%) | **2%** (Amanda, 2026-11-16) | Strictest of the three; on the corrected v4 run, 2% vs 10% changes supported detections by 0.7% (7,270 vs 7,324). |
| 19 | Plasmid hits | count as their host species vs ignore | **Ignore any hit whose title says "plasmid"** (Amanda 2026-11-17) | NCBI files plasmid sequences under the host's taxid, and plasmids move between species; a confirmed C. difficile read's best hit was a C. difficile plasmid. NCBI taxonomy cannot flag these, so this uses the sequence title. |
| 17 | Hits that are not organisms | count them vs ignore them like uncultured | **Ignore** vectors, plasmids, phages and hits filed under "unclassified Bacteria" (named just "bacterium") (Amanda, 2026-11-16) | They are not a competing organism: C. difficile reads tied a phage (Caudoviricetes) on 81% and an "unidentified plasmid" on 75% of confirmed reads; E. faecium reads tied cloning vectors on 85%. Defined by NCBI taxonomy (not under "cellular organisms"), not a name list. |

## B. Still open

| # | Decision point | Options | Status |
|---|---|---|---|
| 20 | 10% rule (advisor): target must be >= 10% of a read's kept named hits (up to 100), else EXCLUDED_NONSPECIFIC | on / off | Built as a switch (off by default); the rerun reports both. Counted over named organism hits only, after dropping uncultured, non-organism and plasmid hits, which limits the bias from species with many genomes in nt. Amanda unsure; decide on the numbers. |
| 16 | Minimum target alignment length | none / 100 bp (reads are 150 bp) | Off by default; the rerun reports 0 and 100 bp. Advisor not keen. A read that would be CONFIRMED on a shorter alignment becomes EXCLUDED_SHORT (neither for nor against). Motivation: B. mallei confirmed reads align over a median 55 bp (all under 100 bp), mostly 16S/conserved fragments shared with unrelated genera. Trade-off at 100 bp: B. mallei 83 -> 0 supported, but V. cholerae 104 -> 52. |
| 15 | Denominator of the validated fraction | CONFIRMED / all reads with a hit vs CONFIRMED / (CONFIRMED + MISASSIGNED) | Both are computed. Only the second lets EXCLUDED / EXCLUDED_COMPLEX reads stop counting against the taxon. |

### Considered and not adopted
| # | Decision point | Why not |
|---|---|---|
| 18 | Target must beat every species from another genus (a cross-genus tie = uninformative) | Too stringent for real taxa: 61% of Salmonella confirmed reads tie E. coli / Klebsiella / Enterobacter (shared Enterobacterales DNA); with #16 Salmonella would fall from 218 to 81 supported detections. Adds little beyond #16 for B. mallei / C. perfringens / V. cholerae. |

## C. Quality checks and what they caught

Each check below found something; each is fixed and covered by tests (40 tests pass).

| Check | Finding | Fix |
|---|---|---|
| Which organism wins the MISASSIGNED reads (`tools/misassigned_winners.py`) | BLAST could not see NCBI's taxonomy file, so every hit was named from its free-text title. 44% of misassigned reads (top 30 taxa) were the correct species under another name (Raoultella/Klebsiella, "[Mycobacterium] chelonae", "MAG:" titles). | BLAST now finds taxdb and also reports taxids (#7); synonym table (#8). |
| Synonym table audit against all NCBI species | Old synonyms that are now other species' names would merge distinct species. | Filter rules (#8). |
| Full v4 rerun, 250 samples (2026-11-11) | All 250 complete; every hit has an NCBI name and taxid. Misassigned reads at 2%: 99,096 (v3) -> 56,837. Raoultella validated fraction 0.21 / 0.33 -> 0.99 / 0.97. | Lab folder swapped to v4; v2 archived. |
| BLAST output of reads counted CONFIRMED, 6 taxa (`make_confirmed_reads_figures.py`) | B. mallei confirmed reads are short (median 55 bp) and tie unrelated genera; C. difficile / E. faecium reads tie phages, plasmids and vectors. | #16 and #17. |
| Pilot rerun, sample GKWWBAA00014 (39 detections, 19,605 reads) | Every hit now has an NCBI name and taxid. Old BLAST + synonym table and new BLAST + taxids give **identical CONFIRMED counts** at 2/5/10% (19,348 / 19,418 / 19,463). MISASSIGNED at 2%: 1,832 before the fix, 226 after (+15 EXCLUDED_COMPLEX). | Two independent routes agree, so the fix is confirmed. |

## D. Next steps
1. Re-score the v4 BLAST output (no new BLAST): 2%, %id >= 90, #17 and #19 on; report the 10% rule and the 100 bp rule on vs off.
2. Decide #20 and #16 from those numbers; Hannah confirms #14 and picks #15.
3. Old vs new per taxon on the final rule, then wire it into the Time-Series filter (#13) with the minimum-confirmed-reads option (#12).
