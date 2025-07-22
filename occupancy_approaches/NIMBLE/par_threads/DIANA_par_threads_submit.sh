#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=4
#SBATCH --job-name=dylcar_explore_occ_run_DIANA_NIMBLE_BUTTERFLIES_par_threads_22_07_2025
#SBATCH --output=slurm_%a.out
#SBATCH --time=24:00:00
#SBATCH --mem=30000
#SBATCH --error=%a.err
#SBATCH --account=ceh_generic
#SBATCH --partition=standard
#SBATCH --qos=high

module load jasr

Rscript --vanilla DIANA_par_threads.r