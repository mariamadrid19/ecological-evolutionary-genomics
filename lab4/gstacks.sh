#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name gstacks
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=10:00:00
#SBATCH -o gstacks.%j.out
#SBATCH -A lp_edu_eeg_2024

module load Stacks/2.5-foss-2018a

gstacks -I alignments/ -M info/popmap.tsv -O stacks.ref/ -t 8
