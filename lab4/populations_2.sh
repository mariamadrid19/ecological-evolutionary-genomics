#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name populations2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=10:00:00
#SBATCH -o populations2.%j.out
#SBATCH -A lp_edu_eeg_2024

module load Stacks/2.5-foss-2018a

populations -P stacks.ref -M info/popmap.tsv -O pop_no_filt --threads 8 --hwe --fstats --vcf --genepop
