#!/bin/bash

#SBATCH -p short-serial

#SBATCH --ntasks=1

#SBATCH --mem=30000

#SBATCH --time=23:59:59

#SBATCH --job-name=ANTS

#SBATCH --array=1-39

#SBATCH --output=Outputs/Ants-%j-%a.txt

module add jasr

Rscript nimble_tests/nimble_ants.R