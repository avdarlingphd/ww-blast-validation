#!/bin/bash
# Run the Snakemake CONTROLLER as its own batch job, so the full run survives SSH
# disconnects / closing your laptop. The controller is lightweight -- it just submits
# and monitors the blast_pair jobs -- but it must stay alive for the whole run, so give
# it generous walltime. The actual compute runs in the per-rule jobs it submits.
#
# Submit from the repo root:
#     sbatch run_snakemake.sbatch
# Watch it:
#     squeue -u $USER            # ww-blast-ctl = this controller; the rest are blast jobs
#     tail -f logs/controller_*.log
# Stop the whole run:
#     scancel <this job id>      # cancelling the controller stops new submissions
#
#SBATCH --job-name=ww-blast-ctl
#SBATCH --partition=hsph,sapphire,shared
#SBATCH --account=hhealy_lab
#SBATCH --time=3-00:00:00          # must outlast the whole run; exits early when done
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/controller_%j.log
set -euo pipefail

ENV=/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/ww-blast-validation
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

# Call the env's snakemake directly (no conda activation needed); it uses the SLURM
# executor profile to submit the per-rule jobs. `sbatch` is on PATH inside a job.
"$ENV/bin/snakemake" -p --profile profiles/slurm
