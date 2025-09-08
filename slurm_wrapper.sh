#!/bin/bash

#SBATCH -A naiss2025-22-471
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 0-12:00:00
#SBATCH -J nf-asm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# This line is specific to the dardel cluster, loads the dependencies into path
ml PDC/24.11 apptainer nextflow

profile=slurm

# run the workflow
nextflow run run_hifiasm.nf -profile ${profile} --input_reads input_reads.csv

## simply add -resume to resume a failed run
#nextflow run run_hifiasm.nf -profile ${profile} --input_reads input_reads.csv -resume