#!/bin/bash
#
#SBATCH --array=0-38
#SBATCH --cpus-per-task=6
#SBATCH --job-name=run_parallel_ants
#SBATCH --output=slurm_%a.out
#SBATCH --time=24:00:00
#SBATCH --mem=40000
#SBATCH --partition=par-single
#SBATCH --error=%a.err
/apps/jasmin/jaspy/miniforge_envs/jasr4.3/mf3-23.11.0-0/envs/jasr4.3-mf3-23.11.0-0-v20240815/lib/R/bin/Rscript --vanilla DIANA_par_species.r
