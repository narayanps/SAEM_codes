#!/bin/bash
#SBATCH --job-name=ecog
#SBATCH --time=48:00:00
#SBATCH --array=11-20,31-40,61-70,101-110,121-130,131-140,141-150
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH -o ecog_out_%A_%a.txt
#SBATCH -e ecog_err_%A_%a.txt
module load teflon
module load matlab
module load CUDA
srun matlab_multithread -nosplash -r "generate_ecog_simulations_all_subs($SLURM_ARRAY_TASK_ID); exit(0)"
