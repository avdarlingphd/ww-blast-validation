# ww-blast-validation: decisions summary (for review with PI)

## BLAST false-positive thresholding rules
- Per candidate (sample x taxid): extract the taxid's reads, subsample to <=500 read pairs (seqtk), blastn vs NCBI nt (max_target_seqs 100).
- Each read is assigned to its single best hit by bitscore; the hit organism is compared to the expected organism.
- TRUE_POSITIVE: >=80% of reads match the expected organism AND expected ranks in the top-2 best-hit organisms AND >=1 matching read AND >=2 reads in the set.
- FALSE_POSITIVE: <10% of reads match.
- UNCULTURED_DOMINANT: >50% of reads best-match uncultured / environmental / metagenome sequences.
- UNCERTAIN: anything in between; NO_DATA: no reads; BLAST_NOT_RUN: not processed.
- Species-specific matching: a species detection needs genus + species-epithet agreement (a genus-level or sibling-species hit does NOT confirm it); genus-level candidates match at genus.

## Directory / data rules
- Code (scripts, Snakefile, config, small resource tables) lives in the GitHub repo (avdarlingphd/ww-blast-validation).
- Large data and outputs (extracted FASTAs, BLAST tables, reports) live on holylabs and are NOT in git (.gitignore excludes them).
- Outputs go to holylabs (/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/blast_validation_v2), NOT netscratch (netscratch is purged).
- Single source of truth = the GitHub repo; old loose copies (Linux Scripts/, 4_blastn_validation/) are deprecated v1 and should not be run.

## Version control / reproducibility
- The cluster runs a git CLONE of the repo; run `git pull` before each run; never edit the cluster copy without committing and pushing.
- Flow: edit + commit + push from laptop -> GitHub (source of truth) -> git pull on cluster -> run.
- conda environment.yaml rebuilds the exact tool versions anywhere (also fixes "snakemake missing on cluster").
- Pin the nt database to a dated snapshot (not "latest") so re-runs are reproducible.
- Stamp each run's output with the git commit hash + the config used.
- Tag pipeline versions (e.g., git tag v2) so a result maps to an exact code state.
