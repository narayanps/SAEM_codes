#!/bin/bash
#SBATCH --job-name=type1_bio
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH -p batch
#SBATCH --array=1-4
#SBATCH -o mcmv_out_%A_%a.txt
#SBATCH -e mcmv_err_%A_%a.txt
module load teflon
module load matlab
srun matlab -nosplash -r "eval_performance_connectivity_v2(1,$SLURM_ARRAY_TASK_ID,'mcmv'); exit(0)"
