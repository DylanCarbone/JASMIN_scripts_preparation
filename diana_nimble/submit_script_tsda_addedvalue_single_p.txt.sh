#!/bin/bash

#SBATCH -p par-single

#SBATCH --cpus-per-task=3

#SBATCH --mem=50000

#SBATCH --time=47:59:59

#SBATCH --job-name=TSDA_single_p

#SBATCH --array=1-46

#SBATCH --output=Outputs/TDSA_Single-p_%j-%a.txt

module add jasr

Rscript added_value/03_hpc_1_singlestream_p.R