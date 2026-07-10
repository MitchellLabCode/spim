#!/bin/bash

#SBATCH --job-name=deconGCaMP      		
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --account=pi-npmitchell
#SBATCH --mem=0G
#SBATCH --cpus-per-task=64
#SBATCH --array=0-399 # range of timepoints for job array


bash deconvolve_timepoints.sh $SLURM_ARRAY_TASK_ID
