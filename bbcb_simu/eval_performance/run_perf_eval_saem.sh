#!/bin/bash
#SBATCH --job-name=type1_bio
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH -p batch
#SBATCH --array=1-4
#SBATCH -o saem_out_%A_%a.txt
#SBATCH -e saem_err_%A_%a.txt
module load teflon
module load matlab
srun matlab -nosplash -r "eval_performance_connectivity_v2(1,$SLURM_ARRAY_TASK_ID,'saem'); exit(0)"
