#!/bin/bash

#SBATCH --job-name=decontest	
#SBATCH --partition=amd
#SBATCH --time=36:00:00
#SBATCH --mem=0G
#SBATCH --cpus-per-task=64
#SBATCH --account=pi-npmitchell
#SBATCH --array=0-149 # range of timepoints for job array
#SBATCH --mail-type=END
#SBATCH --mail-user=wjsh@rcc.uchicago.edu


bash deconvolve_timepoints.sh 
