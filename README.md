# ww-blast-validation

A Snakemake workflow that BLAST-validates candidate taxonomic classifications (e.g. Kraken2 or Bowtie2
calls) in metagenomic samples and filters out classifier false positives. It was built for
hospital-wastewater pathogen surveillance, but the logic is general: give it a list of
`(sample, taxid)` candidate detections plus the reads behind them, and it reports which calls a
rigorous BLAST re-check actually supports.

For each **(sample, taxid)** candidate detection, the pipeline pulls the reads the classifier assigned
to that taxid, BLASTs them against a nucleotide database (NCBI `nt`), and labels the call:

| Classification        | Meaning |
|-----------------------|---------|
| `TRUE_POSITIVE`       | Reads' top BLAST hits match the expected organism at high rate (`avg_pct_match ≥ threshold`, expected organism ranks in the top 2). |
| `FALSE_POSITIVE`      | Reads overwhelmingly hit something else (`avg_pct_match < 10%`). |
| `UNCERTAIN`           | Mixed / intermediate evidence. |
| `UNCULTURED_DOMINANT` | Top hits are dominated by uncultured / environmental / metagenome entries (`> 50%`). |
| `NO_DATA`             | No reads were extracted for this taxid, so there was nothing to BLAST (empty/missing FASTA). |

This removes taxonomic-classifier (Kraken2, or later Bowtie2) false positives from your detection
tables before downstream analysis.

## Pipeline

```
extract_sample  ->  blast_sample  ->  classify
```

1. **`extract_sample`** (method-specific, **batched per sample**) — regenerate the per-(sample, taxid)
   `R1`/`R2` FASTAs of reads the classifier assigned to each taxid. **One job per sample** reads the
   sample's ~5 GB classifier output **once** and both FASTQs **once**, and writes the FASTAs for *all*
   of that sample's taxids in a single pass. A per-sample sentinel (`flags/{sample}.extracted`) tells
   the DAG the group is done.
   - *kraken variant* (`workflow/scripts/extract_kraken_sample.py`): parses `{sample}.kreport` into a
     taxonomy tree and, for each target taxid, collects its whole subtree (**include-children**, so
     strain/subspecies reads under the taxid are captured — faithful to KrakenTools
     `extract_kraken_reads.py --include-children`). It then streams `{sample}.kraken` once to map
     reads → taxids, and streams each FASTQ once to write per-taxid FASTAs. A read whose taxid falls
     under several targets (nested genus + species) is written to each, exactly as running the tool
     once per taxid would.
   - Zero-read taxids still get **empty** FASTAs, so the DAG never breaks; they flow through to
     `NO_DATA`.
   - *Why batched:* the previous per-(sample, taxid) rule reloaded the ~5 GB `.kraken` once per taxid
     (~8,600 loads). Batching to one pass per sample (~250 loads) cuts the `.kraken` and FASTQ read
     I/O by ~34× (the mean taxids/sample) — roughly 43 TB → 1.3 TB of `.kraken` reads over a full run.
2. **`blast_sample`** (**one job per sample**) — subsample each taxid's reads to `max_reads`
   (`seqtk sample -s 42`), **tag every read header with its taxid** (`<taxid>|<R1|R2>|<id>`),
   concatenate all of the sample's candidate reads, and run **one** `megablast` vs `nt` (`-task megablast`,
   `-max_target_seqs 100`, bitscore in the outfmt). This loads the ~1 TB nt DB **once per sample (~250×)**
   instead of once per (sample, taxid) pair (~8,600×) — the DB load dominates runtime, so this is the
   big speedup. Output: one `<sample>_blast_results.tsv`. Empty/missing input → empty result → `NO_DATA`.
3. **`classify`** — split each per-sample TSV back out by the taxid tag and classify every
   (sample, taxid) into `results/<detection_source>/blast_false_positive_report.tsv`.

## Kraken vs Bowtie — a config switch, not a fork

`detection_source` (`kraken` | `bowtie`) is the only thing that differs between candidate sets:

- It selects the method-specific **`extract`** rule (how the per-taxid reads are pulled).
- It routes the **entire output tree** into `results/<detection_source>/` so the two candidate sets
  never collide (`results/kraken/…`, `results/bowtie/…`).

Everything downstream of extraction — `blast_sample`, `classify`, the classifier logic — is
source-agnostic. Adding Bowtie later means writing one `extract_sample` variant that emits the same
output paths + per-sample sentinel; no other rule changes. (The `bowtie` branch is currently a stub.)

## Classification rules

All thresholds live in `config/config.yaml`.

- **Top hit per read, by bitscore.** Each read is judged by its single best BLAST hit (max `bitscore`),
  not by any hit in the returned set.
- **Species-level matching.** A species-level candidate requires the hit to agree at BOTH genus and
  species epithet; a genus-only or sibling-species hit does not confirm it. Genus-level candidates match
  at the genus.
- **`TRUE_POSITIVE`** needs `avg_pct_match ≥ threshold` (default 80), the expected organism in the top 2
  best-hit organisms, and more than `min_match_reads` matching reads per set (default `3`, i.e. > 2).
  A single-read call is never positive.
- **`FALSE_POSITIVE`** is `avg_pct_match < 10%`; **`UNCULTURED_DOMINANT`** is `> 50%` uncultured /
  environmental hits; anything in between is **`UNCERTAIN`**.
- **`-max_target_seqs 100`** so the true organism is not crowded out of the hit list.

See `CHANGELOG.md` for version history.

## Inputs

Set all paths and parameters in `config/config.yaml`. Per sample the pipeline needs:

- the classifier's per-read output + report (Kraken2 `.kraken` / `.kreport`);
- the QC'd paired FASTQs (`{sample}_R1.fastq.gz` / `_R2.fastq.gz`);
- a BLAST nucleotide database (`db:` — point at your own `nt`);

plus two small, in-repo tables: the candidate `(sample, taxid)` list and a `taxid → organism name` map
(`resources/`). Concrete values for the Healy-lab / FASRC deployment are in the example table below.

## Layout

```
config/config.yaml            paths, params, detection_source switch
workflow/Snakefile            rules: extract_sample (kraken) -> blast_sample -> classify
workflow/scripts/
  extract_kraken_sample.py    kraken-variant per-sample batched read extraction (single pass)
  blast_sample.sh             subsample + tag by taxid + one megablast per sample
  false_positive_detection.py v2.1 classifier (splits per-sample TSV by taxid tag)
profiles/slurm/config.yaml    Snakemake SLURM executor profile (FASRC / cannon)
resources/                    candidate (sample, taxid) list + taxid->name map
```

## Running (on the cluster)

Dry run / validate the DAG (no executor or plugin needed):

```bash
snakemake -n
```

Full run via the SLURM profile:

```bash
snakemake -p --profile profiles/slurm
```

The SLURM profile targets FASRC partitions and needs the Snakemake SLURM executor plugin
(`snakemake-executor-plugin-slurm`) installed alongside Snakemake for a real run.

## Testing

Unit tests for the classification logic (species-level matching, the TRUE/FALSE/UNCERTAIN/
UNCULTURED_DOMINANT thresholds, the `>2` matching-reads rule) live in `tests/` and need no cluster:

```bash
pytest tests/                   # with pytest installed (it's in the env)
python tests/test_classify.py   # no pytest needed (self-running)
```

## Example configuration (Healy lab / FASRC)

The concrete paths for our deployment — substitute your own in `config/config.yaml`. Code lives in git;
large data and outputs live on **holylabs** (persistent) and are gitignored.

| What | Where |
|------|-------|
| Classifier inputs (`.kraken` / `.kreport`) | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/kraken_out/kraken_output_ct0_5_min_hit_3/` |
| QC'd paired FASTQs | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip_fastqs/` |
| **Outputs** (FASTAs, BLAST tables, report) | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_blast_validation_v2/` — a **writable sibling** of the read-only `ynhh_ww_rpip_2024/` project dir (keeps raw inputs locked) |
| BLAST `nt` database | **configurable** (`db:`) — point at your local nt; prefer a dated build over `latest`. The build actually used is recorded per run in `PROVENANCE.txt`. |
| Conda environment | `workflow/envs/environment.yaml` |
| Candidate list + names | `resources/` (in git) |

Outputs go to **holylabs, not netscratch** — netscratch is purge-prone and the earlier working tree was
already lost to scratch retention. Nothing under these data paths is committed to git.

## Version control & reproducibility

**The GitHub repo (`avdarlingphd/ww-blast-validation`) is the single source of truth.** The cluster runs
a *clone* of it, never a hand-copied folder — that is what prevents running a stale script.

- First time on the cluster: `git clone https://github.com/avdarlingphd/ww-blast-validation.git`
- Before every run on the cluster: `git pull` (so the clone matches GitHub), then `git status` to
  confirm `working tree clean` / `up to date with origin/main`.
- Edit → commit → push from wherever you work; never edit the cluster copy without committing + pushing
  it back, or the two drift apart.
- Rebuild the environment anywhere: `conda env create -f workflow/envs/environment.yaml`.
- Keep the `nt` database **configurable** (`db:`) — point it at your own local nt rather than relying on
  the FASRC path, and prefer a dated build over `latest`.
- **Automatic provenance.** Each run writes two files next to the report in `results_dir/<source>/`:
  - `PROVENANCE.txt` — the git commit (+ `git describe` and any uncommitted-file list) and the nt database
    build from `blastdbcmd -info`;
  - `config.snapshot.yaml` — a copy of the config as run.

  So every report traces to the exact code, settings, and reference database that produced it — without
  hard-coding any machine-specific path.
- Tag pipeline versions (`git tag v2 && git push --tags`) so a result maps to a fixed code state.

> **Cost warning.** A full run is **250 `extract_sample` jobs + 250 `blast_sample` jobs** against
> `nt` + 1 `classify` (~501 jobs total). Both the `.kraken` reads *and* the ~1 TB nt DB load happen
> once per sample (~250×), not once per taxid (~8,600×) — that batching is what makes the run feasible.
> The `megablast` jobs are still the expensive part; launch the full compute via `sbatch run_snakemake.sbatch`
> (see below) and not without intent.
