#!/usr/bin/env bash
# #####################################################################################################################
# This script calls a Fiji macro deconvolve_timepoints and passes the timepoints to deconvolve.
# Specify start and end timepoints.
#
#
# Example usage
# -------------
# $ bash deconvolve_timepoints.sh 
#
# NPMitchell 2019
#
# Edited and adapted to Midway3 at UChicago RCC, Chris Anto 2024, canto@uchicago.edu
# #####################################################################################################################
source ~/.bashrc
module load fiji/154

echo "executing deconvolution"
xvfb-run fiji -Dimage.updater.disableAutocheck=true --ij2 --run batch_singleTPsingleCh.ijm

# BACK_PID=$!
# wait $BACK_PID
echo "done with deconvolution ch1"

sleep 20s



