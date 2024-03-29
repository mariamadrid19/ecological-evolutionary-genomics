#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name populations
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=10:00:00
#SBATCH -o populations.%j.out
#SBATCH -A lp_edu_eeg_2024

module load Stacks/2.5-foss-2018a

populations -P stacks.ref -M info/popmap.tsv -r 0.8 --min-maf 0.05 --max-obs-het 0.7 -O pop_r0.8_maf0.05_het0.7 --threads 4 --hwe --fstats --vcf --genepop
