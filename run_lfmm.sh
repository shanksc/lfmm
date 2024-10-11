#!/bin/bash
#SBATCH --job-name=lfmm_eur_no_cooper
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --mem=512GB
#SBATCH --cpus-per-task=128
#SBATCH --output=lfmm_eur_no_cooper_k6.log
#SBATCH --time=08:00:00

pwd; hostname; date

snakemake --cores 128 -s prep_lfmm.smk --configfile hg38.yaml --rerun-incomplete

date
