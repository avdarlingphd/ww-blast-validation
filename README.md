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
extract  ->  blast_pair  ->  classify
```

1. **`extract`** (method-specific) — regenerate the per-(sample, taxid) `R1`/`R2` FASTAs of reads the
   classifier assigned to that taxid.
   - *kraken variant:* `KrakenTools/extract_kraken_reads.py -k {sample}.kraken -r {sample}.kreport
     -s R1 -s2 R2 -t {taxid} --include-children --fastq-output`, then `seqtk seq -a` to FASTA.
     `--include-children` walks the `.kreport` taxonomy tree so strain/subspecies reads under the
     taxid are captured.
   - Zero-read taxids still get **empty** FASTAs, so the DAG never breaks; they flow through to
     `NO_DATA`.
2. **`blast_pair`** — per pair: subsample each read set to `max_reads` (`seqtk sample -s 42`), then
   `blastn` vs `nt` with `-max_target_seqs 100` and bitscore in the outfmt.
3. **`classify`** — aggregate every pair's BLAST results into
   `results/<detection_source>/blast_false_positive_report.tsv`.

## Kraken vs Bowtie — a config switch, not a fork

`detection_source` (`kraken` | `bowtie`) is the only thing that differs between candidate sets:

- It selects the method-specific **`extract`** rule (how the per-taxid reads are pulled).
- It routes the **entire output tree** into `results/<detection_source>/` so the two candidate sets
  never collide (`results/kraken/…`, `results/bowtie/…`).

Everything downstream of extraction — `blast_pair`, `classify`, the classifier logic — is
source-agnostic. Adding Bowtie later means writing one `extract` variant that emits the same output
paths; no other rule changes. (The `bowtie` branch is currently a stub.)

## v2 refinements (vs the original per-job shell scripts)

- **Top-hit-by-bitscore.** Each read is judged by its single best hit (max `bitscore`), not "any hit
  in the returned set". `bitscore` is now in the BLAST outfmt.
- **Rescue of high-%match / low-read pairs.** A pair with high `pct_match` but few matching reads is
  now `TRUE_POSITIVE` (a handful of reads that all hit the right organism is confident evidence),
  instead of `UNCERTAIN`. Controlled by `min_match_reads` (default `1`; was `>2`).
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
workflow/Snakefile            rules: extract (kraken) -> blast_pair -> classify
workflow/scripts/
  extract_kraken.sh           kraken-variant read extraction (KrakenTools + seqtk)
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

> **Cost warning.** A full run is ~8,600 extraction jobs + ~8,600 `blastn` jobs against `nt`. The
> kraken variant loads the ~5 GB `.kraken` file once per taxid job; batching per sample would cut that
> I/O. Do not launch the full compute without intent.
