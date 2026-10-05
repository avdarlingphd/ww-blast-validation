# Methods: BLAST false-positive validation (for the manuscript)

Each candidate detection (one sample x taxid from the Kraken2 candidate list) was validated against
the NCBI nt database. For each detection we extracted the reads the classifier assigned to that taxid,
subsampled to at most 500 read pairs (seqtk), and aligned each read to nt with blastn (max_target_seqs
100). Each read was assigned to its single best hit by bitscore, and that hit's organism compared to the
expected organism. For each read set we computed: percent-match (fraction of reads whose best hit is the
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

Parameters: max_target_seqs = 100 (so the true organism is not crowded out of the returned hits);
max_reads = 500 (subsample cap per read set); threshold = 80% (TRUE_POSITIVE); min_match_reads = 3
(more than two matching reads required per set).
