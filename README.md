# ww-blast-validation

A Snakemake workflow that BLAST-validates candidate pathogen detections in hospital wastewater
metagenomes and filters out taxonomic-classifier false positives.

For each **(sample, taxid)** candidate detection, the pipeline pulls the reads the classifier assigned
to that taxid, BLASTs them against NCBI `nt`, and labels the call:

| Classification        | Meaning |
|-----------------------|---------|
| `TRUE_POSITIVE`       | Reads' top BLAST hits match the expected organism at high rate (`avg_pct_match ≥ threshold`, expected organism ranks in the top 2). |
| `FALSE_POSITIVE`      | Reads overwhelmingly hit something else (`avg_pct_match < 10%`). |
| `UNCERTAIN`           | Mixed / intermediate evidence. |
| `UNCULTURED_DOMINANT` | Top hits are dominated by uncultured / environmental / metagenome entries (`> 50%`). |
| `NO_DATA`             | No reads were extracted for this taxid, so there was nothing to BLAST (empty/missing FASTA). |

This is how Kraken2 (or, later, Bowtie2) false positives are removed from the wastewater pathogen
tables.

## Pipeline

```
extract_sample  ->  blast_pair  ->  classify
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
2. **`blast_pair`** — per pair: subsample each read set to `max_reads` (`seqtk sample -s 42`), then
   `blastn` vs `nt` with `-max_target_seqs 100` and bitscore in the outfmt. Depends on the sample's
   extraction sentinel; reads the per-pair FASTA (empty/missing → empty result → `NO_DATA`).
3. **`classify`** — aggregate every pair's BLAST results into
   `results/<detection_source>/blast_false_positive_report.tsv`.

## Kraken vs Bowtie — a config switch, not a fork

`detection_source` (`kraken` | `bowtie`) is the only thing that differs between candidate sets:

- It selects the method-specific **`extract`** rule (how the per-taxid reads are pulled).
- It routes the **entire output tree** into `results/<detection_source>/` so the two candidate sets
  never collide (`results/kraken/…`, `results/bowtie/…`).

Everything downstream of extraction — `blast_pair`, `classify`, the classifier logic — is
source-agnostic. Adding Bowtie later means writing one `extract_sample` variant that emits the same
output paths + per-sample sentinel; no other rule changes. (The `bowtie` branch is currently a stub.)

## v2 refinements (vs the original per-job shell scripts)

- **Top-hit-by-bitscore.** Each read is judged by its single best hit (max `bitscore`), not "any hit
  in the returned set". `bitscore` is now in the BLAST outfmt.
- **Minimum matching reads = `>2`.** `TRUE_POSITIVE` requires more than 2 matching reads per set
  (`min_match_reads: 3`). The v2 single-read "rescue" (`min_match_reads: 1`) was reverted in v2.1 —
  a handful of matching reads is not enough on its own to confirm a detection.
- **`-max_target_seqs 100`** (was `20`) so the true organism is not crowded out of the hit list.
- **Extraction is now part of the workflow.** The old approach BLASTed pre-extracted FASTAs sitting on
  scratch; those were purged by scratch retention, so `extract` regenerates them from the surviving
  Kraken2 outputs + FASTQs.

## Inputs and data provenance

Persistent inputs live on **holylabs** (the netscratch working copies were purged by scratch
retention):

- Kraken2 per-read output + report: `…/ynhh_ww_rpip_2024/kraken_out/kraken_output_ct0_5_min_hit_3/{sample}.kraken` / `.kreport`
- QC'd/dehosted paired FASTQs: `…/ynhh_ww_rpip_2024/Ginkgo_rpip_fastqs/{sample}_R1.fastq.gz` / `_R2.fastq.gz`
- BLAST `nt`: `/n/netscratch/informatics/Everyone/external_repos/blast/nt/latest/nt`
- Candidate list + names: `resources/sample_taxid_list_all_pathogens.tsv`, `resources/taxid_to_name_all_pathogens.csv`
  (8,622 pairs across 250 samples / 142 pathogen taxids).

All paths and parameters are in `config/config.yaml`.

## Layout

```
config/config.yaml            paths, params, detection_source switch
workflow/Snakefile            rules: extract_sample (kraken) -> blast_pair -> classify
workflow/scripts/
  extract_kraken_sample.py    kraken-variant per-sample batched read extraction (single pass)
  blast_pair.sh               subsample + blastn for one pair
  false_positive_detection.py v2 classifier
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

## Data & compute locations

Code lives in git; large data and outputs live on **holylabs** (persistent) and are gitignored.

| What | Where |
|------|-------|
| Classifier inputs (`.kraken` / `.kreport`) | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/kraken_out/kraken_output_ct0_5_min_hit_3/` |
| QC'd paired FASTQs | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip_fastqs/` |
| **Outputs** (FASTAs, BLAST tables, report) | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/blast_validation_v2/` |
| BLAST `nt` database | pin to a **dated** nt snapshot (see config `db:` TODO), not `.../nt/latest/nt` |
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
- Pin the `nt` database to a dated snapshot (config `db:`), not `latest`.
- Stamp each run's output with the commit it ran: `git rev-parse HEAD > <results_dir>/COMMIT.txt`
  and keep a copy of `config/config.yaml` alongside the report, so any result traces to exact code + settings.
- Tag pipeline versions (`git tag v2 && git push --tags`) so a result maps to a fixed code state.

### Deprecated copies

Older loose copies of `false_positive_detection.py` — under `.../Linux Scripts/` and
`.../Biofilm-Project/4_blastn_validation/` — are **deprecated v1** and must not be run. They have drifted
from this repo; **this repository is canonical**.

> **Cost warning.** A full run is **250 `extract_sample` jobs + ~8,600 `blast_pair` jobs** against
> `nt` + 1 `classify` (~8,875 jobs total). Extraction now reads each ~5 GB `.kraken` once per sample
> (~250 loads) instead of once per taxid (~8,600 loads). The `blastn` jobs are the expensive part;
> do not launch the full compute without intent.
