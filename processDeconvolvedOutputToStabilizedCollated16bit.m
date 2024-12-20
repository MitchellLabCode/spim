%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
cd /mnt/crunch/Chris/5-28_experiment_line/Deconvolved_2024-05-28_150855/All_timepoint_channels/
zwidth =1 ;
dataDir = cd ;
%directory name: deconvolved_16bit_pre_pressure_-10_post_pressure_5_TP30-67

%%Save data
%save('/mnt/data/code/haibei/matlab/haibei_scripts/midgutClosure'); %save workspace
%save('/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/deconvolved_16bit_NEW/midgut_tubulogenesis_workspace_variables');
%% ADD PATHS TO THIS ENVIRONMENT ==========================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
%Do not do this. Just open the tubular folder in my folder
% addpath(fileparts(origpath))
%addpath(genpath('/mnt/data/code_static/tubular/'))
%addpath('/mnt/data/code_static/tubular/TexturePatch


addpath(genpath('/mnt/data/code/gptoolbox'))
addpath(genpath('/mnt/data/code/canto/tubular'))
% go back to the data
cd(dataDir)




%% DEFINE NEW MASTER SETTINGS

% Pay close attention to scale. This parameter defines the max intensity of
% the tissue that you care about. 
if ~exist('./masterSettings.mat', 'file')  
    % Metadata about the experiment
    stackResolution = [.195 .195 .195] ;  % resolution in spaceUnits per pixel %every pixel in space is 0.2619 spaceUnit (um)? 
    nChannels = 3 ;             % how many channels is the data (ex 2 for GFP + RFP)
    channelsUsed = [0 1 2] ;          % which channels are used for analysis
    %timePoints = 1:10:101;       % timepoints to include in the analysis
    %The timepoints should be changed. Everything after tp 71 is too late.
    %So we do not include these. 
    % timePoints = 1:10:71;
    % timePoints = 0:1; 
    timePoints = [0:119];
    %timePoints = [0, 20, 60];

    ssfactor = 5 ;              % subsampling factor
    flipy = true ;             % whether the data is stored inverted relative to real position in lab frame
    timeInterval = 2 ;          % physical interval between timepoints
    timeUnits = 'min' ;         % physical unit of time between timepoints
    spaceUnits = [char(956) 'm'] ;     % physical unit of time between timepoints
    scale = [3 0.1 0.7] ;      % scale for conversion to 16 bit, one value for each channel
    %file32Base = 'TP%d_Ch%d_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    file32Base = 'TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'Time_%06d_c%d_stab';  % stabilized (jitter removed) filename
    fn_prestab = 'Time_%06d_c%d.tif'; % pre-stabilization filename
    %set_preilastikaxisorder = 'xyzc' ; % data axis order for subsampled h5 data (ilastik input)
    %set_preilastikaxisorder = 'yzxc' ;
    set_preilastikaxisorder = 'cxyz' ;
    swapZT = 0 ;                % whether to swap the z and t dimensions
    masterSettings = struct('stackResolution', stackResolution, ...
        'scale', scale, ...
        'file32Base', file32Base, ...
        'fn_prestab', fn_prestab, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', spaceUnits, ...
        'fn', fn,...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        'nU', 100, ...  
        'nV', 100); 
    disp('Saving masterSettings to ./masterSettings.mat')
    if exist('./masterSettings.mat', 'file')
        ui = input('This will overwrite the masterSettings. Proceed (Y/n)?', 's') ;
        if ~isempty(ui) && (strcmp(ui(1), 'Y') || strcmp(ui(1), 'y'))
            save('./masterSettings.mat', 'masterSettings')
            loadMaster = false ;
        else
            disp('Loading masterSettings from disk instead of overwriting')
            loadMaster = true ;
        end
    else
        save('./masterSettings.mat', 'masterSettings')
        loadMaster = false ;
    end
else
    loadMaster = true ;
end

if loadMaster
    % LOAD EXISTING MASTER SETTINGS
    disp('Loading masterSettings from ./masterSettings.mat')
    load('./masterSettings.mat', 'masterSettings')
    % Unpack existing master settings
    stackResolution = masterSettings.stackResolution ;
    nChannels = masterSettings.nChannels ;
    channelsUsed = masterSettings.channelsUsed ;
    timePoints = masterSettings.timePoints ;
    ssfactor = masterSettings.ssfactor ;
    % whether the data is stored inverted relative to real position
    flipy = masterSettings.flipy ; 
    timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
    timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints
    spaceUnits = masterSettings.spaceUnits ; % unit of distance of full resolution data pixels ('$\mu$m')
    fn = masterSettings.fn ;
    fn_prestab = masterSettings.fn_prestab ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
    swapZT = masterSettings.swapZT ;
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
    scale = masterSettings.scale ;
    nChannels = 2 ;             % how many channels is the data (ex 2 for GFP + RFP)
    channelsUsed = [0 1 2] ;
    timePoints = [10 70];
end



%dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;

% dir32bit = fullfile(dataDir, 'deconvolved_12iter18iter_every_10_timepoints') ;
% dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
% dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;


%In stabilizeImages, we have already got the MIPs of the channels. We use
%them to remove jitter. 
%%We have run the code above, however we are not sure the value of
%%im_intensity and imref_intensity? 

dir32bit = fullfile(dataDir, '') ;
% dir32bit_tp50 = fullfile(dataDir, 'deconvolved_12iter_ch2_06132024_tp50');
%dir16bit = fullfile(dataDir, 'deconvolved_16bit_NEW') ;
%dir16bit = fullfile(dataDir, 'deconvolved_16bit_CHECK_AXIS_ORDER/');
dir16bit = fullfile(dataDir, 'deconvolved_16bit_0810');
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;
%% END OF EXPERIMENT METADATA =============================================
% =========================================================================
% =========================================================================
addpath(genpath('/mnt/data/code_static/matlab/plotting/'))
addpath(genpath('/mnt/data/code_static/matlab/tiff_handling/'))
%
% Skip if already done
convert32to16bit(timePoints, scale, dir32bit, dir16bit_prestab,...
file32Base, fn_prestab, channelsUsed)



%% -III. make MIPs for 16bit images
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
% % Rename stab to prestab -- hack
% % fns = fullfile('./deconvolved_16bit/Time*stab')
% % for qq=1:length(fns)
% % command = ['mv ' fullfile(fns.folder, fns.name) fullfile(dir16bit, fns.name
% % end

% -III. make MIPs for 16bit images -- optional
% Skip if already done
mipDir = fullfile(dir16bit_prestab, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)

%% -II. stabilize images, based on script stabilizeImagesCorrect.m
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
Options.im_intensity = 1 ; % 0.01 ; %scale for MIP output intensity
Options.imref_intensity = 1 ; % 0.005 ; %scale for reference MIP intensity
%Specific meaning of these two variables (im_intensity and imref_intensity: See the notes in "stabilizeImages.m")
Options.stabChannel = 0 ; % which channel I am stabilizing. Default = 1.

stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
timePoints, timePoints, timePoints(1), ...
mipDir, mipoutdir, mips_stab_check, Options)

%% Stabilize the other channels based on stab_channel withour rewriting mips
Options.overwrite_mips = false ;
Options.stabChannel = 0 ;
Options.overwrite_tiffs = true ;
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
timePoints, timePoints, timePoints(1), ...
mipDir, mipoutdir, mips_stab_check, Options)

Options.overwrite_mips = false ;
Options.overwrite_tiffs = true ;
Options.stabChannel = 2 ;
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
timePoints, timePoints, timePoints(1), ...
mipDir, mipoutdir, mips_stab_check, Options)
%% Collate multiple colors into one TIFF post-stabilization
%collate: collect and combine in proper order
%timePoints=[33:]
fileNameIn = fullfile(dir16bit, 'Time_%06d_c%d_stab.tif') ;
fileNameOut = fullfile(dir16bit, 'Time_%06d.tif') ;
collateColors(fileNameIn, fileNameOut, timePoints, channelsUsed) ; 

%% Once the data is collated into multi-color tiffs, delete the individual
% channel tiffs
