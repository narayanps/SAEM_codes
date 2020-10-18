#!/bin/bash
#SBATCH --job-name=type1_bio
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH -p batch
#SBATCH --array=1-200
#SBATCH -o bio_int_dip_out_%A_%a.txt
#SBATCH -e bio_int_dip_err_%A_%a.txt
module load teflon
module load matlab
srun matlab -nosplash -r "run_bbcb_sim($SLURM_ARRAY_TASK_ID,1,0); exit(0)"
