# Methods: BLAST false-positive validation (for the manuscript)

Each candidate detection (one sample x taxid from the Kraken2 candidate list) was validated against
the NCBI nt database. For each detection we extracted exactly the reads Kraken2 had assigned to that
taxid (its subtree, via --include-children), subsampled each read direction (R1, R2) to at most 500
reads (seqtk sample -s 42; in practice ~86% of candidates had fewer than 500 reads and so used all of
them), and aligned each read to nt with blastn (megablast task, max_target_seqs 100). Each read was
assigned to its single best hit by bitscore, and that hit's organism compared to the expected organism.
Because only the reads already assigned to the target are BLASTed, and the subsample estimates the
fraction of those reads that match (standard error <~2.5% at n = 500), the cap does not cause the target
to be "missed"; low-abundance targets (<500 reads) are validated on all of their reads. For each read set we computed: percent-match (fraction of reads whose best hit is the
expected organism), expected-rank (rank of the expected organism among best-hit organism frequencies),
and percent-uncultured (fraction of reads best-matching uncultured / environmental / metagenome entries).

Averaging the R1 and R2 read sets, a detection was classified as:
- UNCULTURED_DOMINANT if more than 50% of reads best-matched uncultured/environmental/metagenome sequences;
- TRUE_POSITIVE if at least 80% of reads matched the expected organism, the expected organism ranked in
  the top two best-hit organisms, and more than two reads matched the expected organism (min_match_reads = 3;
  a handful of matching reads is not sufficient and is called UNCERTAIN);
- FALSE_POSITIVE if fewer than 10% of reads matched the expected organism;
- UNCERTAIN otherwise; NO_DATA if no reads were recovered; BLAST_NOT_RUN if the pair was not processed.

Organism matching was species-specific: a species-level candidate required agreement of BOTH genus and
species epithet, so that a genus-level hit or a sibling-species hit did not confirm a species-level
detection; genus-level candidates were matched at the genus level.

Hit ranking. For each read, BLAST returns up to max_target_seqs alignments ordered by E-value
(ascending) — the E-value being the number of alignments expected by chance in a database of this size
scoring at least as well, so lower is more significant. We select each read's single best hit by
bitscore rather than by raw percent identity: bitscore accounts for both alignment identity and length
(and is normalized to database size), so a long near-perfect alignment outranks a short perfect one, and
within a single query bitscore is monotonic with E-value. max_target_seqs was set to 100 because the
parameter is a collection heuristic, not a guaranteed top-N by score — BLAST retains the first N targets
it encounters during the search and then sorts them, so a low value can omit the true best hit (Shah et
al. 2019, Bioinformatics 35:1613-1614). A value of 100 ensures the genuine best hit is captured and
provides the competing-organism context used for the expected-rank and top-organisms fields.

Parameters: max_target_seqs = 100 (see Hit ranking); max_reads = 500 (subsample cap per read set);
threshold = 80% (TRUE_POSITIVE); min_match_reads = 3 (more than two matching reads required per set).
The nt database build used for each run was recorded automatically (blastdbcmd -info) alongside the git
commit of the code, so results trace to an exact code + reference-database version.
