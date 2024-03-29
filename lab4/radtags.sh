#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name radtags
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=10:00:00
#SBATCH -o radtags.%j.out
#SBATCH -A lp_edu_eeg_2024

module load Stacks/2.5-foss-2018a

process_radtags -1 /data/leuven/357/vsc35707/lab4/raw/lib1.R1.fastq.gz -2 /data/leuven/357/vsc35707/lab4/raw/lib1.R2.fastq.gz -b /data/leuven/357/vsc35707/lab4/info/barcodes.lib1.tsv --inline-inline -o /data/leuven/357/vsc35707/lab4/cleaned/ --renz-1 sbfI --renz-2 sphI --rescue --clean --quality --adapter_1 CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT --adapter_2 AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter_mm 2
