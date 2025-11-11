%% Make MIPs for either fused or deconvolved images


datdir = '/project/npmitchell/avistrok/bapGAL4_UAShidUASstingerHiRFP/2025-08-04_131451/unpacked/fused/'

% /project/npmitchell/canto/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-29_232230_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/unpa
fns = dir(fullfile(datdir, '*.tif')); 
% get timepoint list from list of timepoints already done
fnPattern = 'TP%0d_Ch%d_Ill0_Ang4,64,124,184,244,304.tif';
% todo: get list from pattern
tpList = [0:76] ;

for i=1:length(tpList)
    fname=fns(i).name;
    tokens=regexp(fname,'TP(\d+)_Ch0','tokens');
    if ~isempty(tokens)
        tp=str2double(tokens{1}{1});
        tpList(end+1)=tp;

    end
end

timePoints=unique(tpList);

addpath('/project/npmitchell/canto/code/matlab/master_pipeline/')
addpath('/project/npmitchell/canto/code/tubular/utility/bfmatlab/')

outdir = fullfile(datdir, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
Options.channels = [0] ;
makeMips(timePoints, datdir, fnPattern, outdir, Options)

%/project/npmitchell/canto/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-29_232230_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/unpa
fns = dir(fullfile(datdir, '*.tif')); 
% get timepoint list from list of timepoints already done
fnPattern = 'TP%0d_Ch%d_Ill0_Ang4,64,124,184,244,304.tif';
% todo: get list from pattern
tpList = [0:76] ;

for i=1:length(tpList)
    fname=fns(i).name;
    tokens=regexp(fname,'TP(\d+)_Ch1','tokens');
    if ~isempty(tokens)
        tp=str2double(tokens{1}{1});
        tpList(end+1)=tp;

    end
end

timePoints=unique(tpList);

addpath('/project/npmitchell/canto/code/matlab/master_pipeline/')
addpath('/project/npmitchell/canto/code/tubular/utility/bfmatlab/')

outdir = fullfile(datdir, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
Options.channels = [1] ;
makeMips(timePoints, datdir, fnPattern, outdir, Options)

% 
% % /project/npmitchell/canto/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-29_232230_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/unpa
% fns = dir(fullfile(datdir, '*.tif')); 
% % get timepoint list from list of timepoints already done
% fnPattern = 'TP%0d_Ch%d_Ill0_Ang36,96,156,216,276,336.tif';
% % todo: get list from pattern
% tpList = [0:52] ;
% 
% for i=1:length(tpList)
%     fname=fns(i).name;
%     tokens=regexp(fname,'TP(\d+)_Ch2','tokens');
%     if ~isempty(tokens)
%         tp=str2double(tokens{1}{1});
%         tpList(end+1)=tp;
% 
%     end
% end

% timePoints=unique(tpList);
% 
% addpath('/project/npmitchell/canto/code/matlab/master_pipeline/')
% addpath('/project/npmitchell/canto/code/tubular/utility/bfmatlab/')
% 
% outdir = fullfile(datdir, 'mips') ;
% Options.overwrite_mips = false ;
% Options.scale = -1 ; % do NOT rescale intensities during intensity projection
% Options.channels = [2] ;
% makeMips(timePoints, datdir, fnPattern, outdir, Options)
