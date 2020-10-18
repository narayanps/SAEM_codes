#!/bin/bash
#SBATCH -p batch
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --array=1-19
#SBATCH -o preproc_out_%A_%a.txt
#SBATCH -e preproc_err_%A_%a.txt
module load anaconda3
module load teflon
srun python run_meg_preprocessing_v2.py $SLURM_ARRAY_TASK_ID
