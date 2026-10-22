# Design decisions & reasoning — BLAST validation classifier

Running record of the methodological choices and the reasoning behind them (for the manuscript, for
advisor discussions, and for continuity). Dated as decisions are made.

## Subsampling: keep the ≤500-reads-per-target cap  (decided 2026-10)
- **Scope:** subsampling is per **(target taxid, sample)** — each pathogen's reads *within* a sample are
  capped at 500 independently (not one cap per sample). seqtk subsamples **whole reads** (never parts of a
  read); every selected read is BLASTed in full.
- **Why keep it:** `pct_match` is a proportion, compared to coarse thresholds (10% / 80%). 500 reads
  estimate that proportion to a standard error of ~±2% — far tighter than the thresholds need. ~86% of
  targets have <500 reads (all reads used); only ~13.5% hit the cap.
- **The capped targets are the abundant, well-identified organisms** (P. aeruginosa, Aeromonas caviae,
  H. parainfluenzae, Morganella morganii, B. fragilis, E. faecium, K. quasipneumoniae, Enterobacter
  cloacae, E. coli, …) and **87% of capped detections are already TRUE_POSITIVE on 500 reads.** More reads
  would re-confirm confident calls at large compute cost and change ~nothing.
- Removing the cap is an *upstream* change → full re-BLAST, big query-size blow-up, timeout risk, ~zero
  changed calls. Not worth it. (Methods wording: "subsampled to a fixed depth sufficient to estimate the
  per-target match proportion to ±~2%.")

## Hit ranking & max_target_seqs = 100
- `-max_target_seqs 100` is **per read** (each ~150 bp read is its own query), returning up to 100
  **database accessions** (redundant — many per organism), not 100 organisms.
- Hits are ordered by **E-value ascending ≈ bitscore descending**. We pick each read's best hit by
  **bitscore** (not raw % identity): bitscore combines identity + alignment length + database size, so a
  long near-perfect alignment outranks a short perfect one.
- 100 is used because `max_target_seqs` is a *collection heuristic* (first-N-found, then sorted; Shah et
  al. 2019, Bioinformatics 35:1613), so a low value can drop the true best hit; 100 also captures the
  tied top cluster.

## Rejected: "≥50% of the top-100 hits must be the target"
- The top-100 are **redundant accessions**, so the percentage is confounded by nt **deposition bias**: an
  under-represented true organism (e.g. *K. variicola*, few nt entries) can be <50% of its own read's
  top-100 while a tie-related over-represented sibling (*K. pneumoniae*, tens of thousands of entries)
  dominates — the rule would wrongly fail a correct call, and inflate over-represented ones. With ties,
  it mostly measures "which organism has more deposited copies." Use distinct-organism ranking instead.

## Adopt: tie-aware matching  (preferred; implement as optional flag)
- **Current:** per read, take the single #1 hit by bitscore; on a tie, sort order decides (effectively
  arbitrary, and biased toward the more-deposited sibling). `pct_match` = fraction of reads whose #1 = target.
- **Tie-aware:** a read counts for the target if the target reaches that read's **maximum** bitscore
  (ties count). De-duplicated by organism, database-bias-free, per-read. Fixes the sibling-tie distortion
  (worked example: 40% → 80% for a *K. variicola* candidate). Reads that match a *different* organism
  strictly better are still excluded, so it is not over-permissive.
- **vs `expected_rank`:** `expected_rank ≤ 2` is a *set-level frequency* check ("is the target among the
  top-2 most-common per-read winners across the read set?"); tie-aware `pct_match` is a *per-read*
  membership check ("on this read, is the target a best hit?"). Different levels; the per-read one is what
  corrects the tie distortion.

## Adopt: treat uncultured best hits as "uninformative"  (liked; implement + test)
- **Current:** a read whose best hit is uncultured/environmental/metagenome counts toward
  `pct_uncultured`; >50% → UNCULTURED_DOMINANT (checked first).
- **New idea:** treat an uncultured *best hit* as **uninformative** — exclude it from the denominator
  (neither match nor non-match) rather than counting it as a non-match. Many "uncultured bacterium clone"
  nt entries *are* the organism under an uncultured label, so counting them as non-matches under-credits a
  real target. Keep UNCULTURED_DOMINANT as a separate flag so genuinely uncultured-dominated calls are
  still surfaced.

## Clarification: a sibling tie is NOT a false positive
- A read tying the target and a complex sibling (*K. variicola* ↔ *K. pneumoniae*) = **species ambiguity**
  (genuine genus/complex, species unresolvable), not a false positive. A false positive is when reads
  best-hit a **different lineage** entirely. The complex-level scoring handles the former.

## Open / to build
- [ ] Implement **tie-aware matching** option; run on the saved per-read TSVs (free, no re-BLAST); report
      how many detections move and which organisms (expect the Klebsiella/Enterobacter/Acinetobacter siblings).
- [ ] Implement **uncultured-uninformative** option; test the same way.
- [ ] **Clade-level capture:** add genus/complex targets, re-extract + re-BLAST, to validate reads Kraken
      parked at a higher clade (rescues NO_DATA/ambiguous species like *Shigella boydii*). Separate re-run.
- [ ] Replace the **provisional** complex crosswalk with the real **RPIP panel v1** and **YNHH clinical**
      complex schemes once available.
- [ ] **Threshold 75 vs 80** sensitivity pass (free re-score of saved TSVs).
