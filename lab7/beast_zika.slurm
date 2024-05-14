#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --time=00:010:00
#SBATCH --account=lp_edu_eeg_2024
#SBATCH --cluster=genius

module load Java
module load BEAGLE
module load Beast/1.10.4-intel-2018a

beast zika.xml
