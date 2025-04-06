%% Make MIPs for either fused or deconvolved images


datdir = './' ;
fns = dir(datdir) ;
% get timepoint list from list of timepoints already done
fnPattern = 'Time_%06d_channel_%a_Angle0,60,120,180,240,300.tif';
% todo: get list from pattern
timePoints = finish_here ;

outdir = fullfile(datdir, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
Options.channels = [0 1 2] ;
makeMips(timePoints, datdir, fn_prestab, outdir, Options)
