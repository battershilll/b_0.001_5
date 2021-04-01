#!/bin/bash -e
#SBATCH --job-name=0.0_1_17.5_visc1
#SBATCH --time=60:00:00 #Walltime (HH:MM:SS)
#SBATCH --mem=20G
#SBATCH --output=out.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24


#SBATCH --profile task 
pwd #Prints working directory

ml use /scale_wlg_persistent/filesets/project/uoa02652/.local/modules/all # Include local modules
ml Basilisk/Matty-gimpi-2018b


./test_robin




