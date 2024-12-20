%% processDeconvolvedOutputToStabilizedCollated16bit.m
% Script for processing output TIFFs from the multi-view deconvolution pipeline
% Things to change each time are: 
% dataDir : path to the output data 
% scale : max intensity values per channel for renormalizing when 
%        downsampling to 16bit
% file32Base : filename pattern for the muvi deconvolution output
% Options.channels : list of channels, ex [0,1,2] for 3 channel imaging 
%        with ch0 ch1 ch2 
% Options.stabChannel : which channel to use for jitter stabilization. 
%        Use a channel where the beads are much brighter than the sample
% referenceTimepoint : which timepoint to use as a reference for others to 
%        be translated to match. 
%        For ex, if the sample is drifting over time, be careful that 
%        subtracting  off this motion doesn't push your ROI outside the 
%        image volume.

%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
dataDir = '/mnt/crunch/Chris/5-28_experiment_line/Deconvolved_2024-05-28_150855/All_timepoint_channels/'
cd(dataDir)
addpath(genpath('/project/npmitchell/code/matlab/'))


%% DEFINE MASTER SETTINGS
scale = [3 0.1 0.7] ;      % scale for conversion to 16 bit, one value for each channel
file32Base = 'TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif'; 
fn = 'Time_%06d_c%d_stab';  % stabilized (jitter removed) filename

dir32bit = fullfile(dataDir, '') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit');
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;
%% END OF EXPERIMENT METADATA =============================================
% =========================================================================

addpath(genpath('/mnt/data/code_static/matlab/plotting/'))
addpath(genpath('/mnt/data/code_static/matlab/tiff_handling/'))

% -V. Convert to 16 bit
convert32to16bit(timePoints, scale, dir32bit, dir16bit_prestab,...
file32Base, fn_prestab, channelsUsed)

%% -IV. make MIPs for 16bit images
% Skip if already done
mipDir = fullfile(dir16bit_prestab, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
Options.channels = [0 1 2] ;
makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)

%% -III. make saggital MIPs for 16bit images -- optional
% Skip if already done
subMipDir = fullfile(dir16bit_prestab, 'substack') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
Options.overlayColors = true ;
Options.width = 5 ;
Options.channels = [0 1 2] ;
Options.overwrite_overlays = false ;
Options.pages = [607:611] ;
Options.dim = 2 ;
makeSubStackMips(timePoints, dir16bit_prestab, fn_prestab, subMipDir, Options)

%% Collate multiple colors into one TIFF pre-stabilization
fileNameIn = fullfile(dir16bit_prestab, fn_prestab) ;
fileNameOut = fullfile(dir16bit_prestab, 'Time_%06d.tif') ;
collateColors(fileNameIn, fileNameOut, timePoints, channelsUsed) ;

% -II. make MIPs for 16bit images -- optional
% Skip if already done
mipDir = fullfile(dir16bit_prestab, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)

%% -I. stabilize images, based on script stabilizeImagesCorrect.m
% There is some jitter from stage motion (or sample motion). We want to
% correct for this. We measure motion of the sample in XYZ and substract
% off that motion. We don't correct for rotations, and that is a limitation
% of this step. But in practice, subtracting off translations is helpful
% and sufficient. This is done by performing phase correlation (ie
% correlation in Fourier space).
% Skip if already done
% name of directory to check the stabilization of mips
mips_stab_check = fullfile(mipDir, 'stab_check') ; %Output layer for the
mipoutdir = fullfile(mipDir, 'mips_stab') ; %Output dir for the stabilized MIPs
% Choose bit depth as typename
typename = 'uint16' ;
% Give file names for I/O
fileNameIn = fullfile(dir16bit_prestab, fn_prestab) ;
fileNameOut = fullfile(dir16bit, [fn '.tif']) ;
rgbName = [fn '.png'] ;
typename = 'uint16' ;
Options = struct();
Options.overwrite_mips = false ;
Options.overwrite_tiffs = true ;
Options.im_intensity = 1 ; 
Options.imref_intensity = 1 ; 
%Specific meaning of these two variables (im_intensity and imref_intensity: See the notes in "stabilizeImages.m")
Options.stabChannel = 0 ; % which channel I am stabilizing. Default = 1.

referenceTimepoint = timePoints(end) ;  % Timepoint that is not affected by stabilization (set as reference)
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
timePoints, timePoints, referenceTimepoint, ...
mipDir, mipoutdir, mips_stab_check, Options)

%% Stabilize the other channels based on stab_channel withour rewriting mips
Options.overwrite_mips = false ;
Options.stabChannel = 1 ;
Options.overwrite_tiffs = true ;
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
    timePoints, timePoints, referenceTimepoint, ...
    mipDir, mipoutdir, mips_stab_check, Options)

Options.overwrite_mips = false ;
Options.overwrite_tiffs = true ;
Options.stabChannel = 2 ;
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
timePoints, timePoints, referenceTimepoint, ...
mipDir, mipoutdir, mips_stab_check, Options)

%% Collate multiple colors into one TIFF post-stabilization
%collate: collect and combine in proper order
fileNameIn = fullfile(dir16bit, 'Time_%06d_c%d_stab.tif') ;
fileNameOut = fullfile(dir16bit, 'Time_%06d.tif') ;
collateColors(fileNameIn, fileNameOut, timePoints, channelsUsed) ; 

%% Once the data is collated into multi-color tiffs, delete the individual
% channel tiffs
