#!/usr/bin/env bash
# #####################################################################################################################
# This script calls a Fiji macro deconvolve_timepoints and passes the timepoints to deconvolve.
# Specify start and end timepoints.
#
#
# Example usage
# -------------
# $ bash deconvolve_timepoints.sh 73
# In this example, 73 is the timepoint to be deconvolved/fused
#
# NPMitchell 2019
#
# Edited and adapted to Midway3 at UChicago RCC, Chris Anto 2024, canto@uchicago.edu
# #####################################################################################################################
source ~/.bashrc
module load fiji

echo "executing deconvolution"
xvfb-run -d /project/npmitchell/fiji-linux64/Fiji.app/ImageJ-linux64 -Dimage.updater.disableAutocheck=true --ij2 --console --run batchtest.bsh $timepointID
BACK_PID=$!
wait $BACK_PID
echo "done with deconvolution"
sleep 200s



