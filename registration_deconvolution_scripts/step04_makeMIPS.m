%% Make MIPs for either fused or deconvolved images


datdir = '/project/npmitchell/canto/HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/2024-06-07_CAAX/unpacked_flipped_rotated_TIFFs_new_06-07-CAAX/deconvolved32bit/' ;
fns = dir(fullfile(datdir, '*.tif')); 
% get timepoint list from list of timepoints already done
fnPattern = 'TP%0d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif';
% todo: get list from pattern
tpList = [] ;

for i=1:length(fns)
    fname=fns(i).name;
    tokens=regexp(fname,'TP(\d+)_Ch0','tokens');
    if ~isempty(tokens)
        tp=str2double(tokens{1}{1});
        tpList(end+1)=tp;

    end
end

timePoints=unique(tpList);



outdir = fullfile(datdir, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
Options.channels = [0] ;
makeMips(timePoints, datdir, fn_prestab, outdir, Options)
