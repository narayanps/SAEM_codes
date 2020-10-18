#!/bin/bash
#SBATCH -p batch
#SBATCH --time=20:00:00      # 4 hours
#SBATCH --mem-per-cpu=32G   # 1G of memory
#SBATCH -o preproc_out_%A_%a.txt
#SBATCH -e preproc_err_%A_%a.txt
module load anaconda3
module load teflon
srun python fetch_data.py
