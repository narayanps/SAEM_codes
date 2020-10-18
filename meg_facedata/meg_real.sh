#!/bin/bash
#SBATCH --job-name=MEG
#SBATCH --time=72:00:00
#SBATCH --mem=20G
#SBATCH --array=101-160
#SBATCH --gres=gpu:1
#SBATCH -o MEG_out_%A_%a.txt
#SBATCH -e MEG_err_%A_%a.txt
module load teflon
module load matlab
module load CUDA
srun matlab_multithread -nosplash -r "run_saem_for_meg($SLURM_ARRAY_TASK_ID, meg_path, save_path); exit(0)"
