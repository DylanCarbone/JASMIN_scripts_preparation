#!/bin/bash

#SBATCH -p par-single

#SBATCH --ntasks=4

#SBATCH --mem=50000

#SBATCH --time=47:59:59

#SBATCH --job-name=NIMBLE_parallel

#SBATCH --array=1

#SBATCH --output=Outputs/Test.txt

module add jasr

Rscript nimble_parallel.R