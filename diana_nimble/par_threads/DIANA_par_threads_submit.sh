#!/bin/bash
#
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --job-name=run_parallel_threads
#SBATCH --output=slurm_%a.out
#SBATCH --time=02:00:00
#SBATCH --mem=20000
#SBATCH --partition=par-single
#SBATCH --error=%a.err

module load jasr

/apps/jasmin/jaspy/miniforge_envs/jasr4.3/mf3-23.11.0-0/envs/jasr4.3-mf3-23.11.0-0-v20240815/lib/R/bin/Rscript --vanilla DIANA_par_threads.r
