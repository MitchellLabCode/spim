%% Setting up ImSAnE
%
% You may need to go into the imsane directory and run setup.m
%

%% Initialize ImSAnE project
%
% We start by clearing the memory and closing all figures.
%

clear all ;
close all ;

%% Experiment Options
%
% Set options for the experiment. The options are listed in sections based
% on what each one is used for.
%

%
% General options for the experiment.
%
% * 'overwrite' , the boolean whether to overwrite previous pullbacks
% * 'visualize' , the boolean whether to show figures
% * 'fnfull'    , the file name convention for fused data
% * 'zDim'      , the long axis of the embryo
%

overwrite = true ;
visualize = true ;
fnfull    = 'Time_%06d_c' ;
zDim      = 1 ;

%
% Options for time and channel detection
%
% * 'firstTime'  , the first time in the data set
% * 'UserTimes'  , the number of sequential timepoints after firstTime
% * 'detecttime' , the boolean whether to detect number of times
% * 'firstChan'  , the first channel in the data set
% * 'UserChan'   , the number of sequential channels after firstChan
% * 'detectchan' , the boolean whether to detect number of channels
%

firstTime  = 5 ;
UserTimes  = 2 ;
detecttime = false ;
firstChan  = 1 ;
UserChan   = 1 ;
detectchan = true ;

%
% Options for the file metadata.
%
% * 'stackResolution' , the stack resolution in microns
% * 'swapZT'          , the option to swap the third and fourth dimension
%

stackResolution = [.2619 .2619 .2619] ;
swapZT          = 1 ;

%
% Options for experiment metadata.
%
% * 'channelsUsed' , the channels used
% * 'channelColor' , the mapping from element in channels used to RGB = 123
% * 'description'  , the description of the data
%

channelsUsed = 1 ;
channelColor = 1 ;
description  = 'Stained embryo' ;

%
% Options for detecting the initial point cloud.
%
% * 'guessChannel'             , the integer specifying which channel(s) to
%                                use for detection
% * 'guessSigma'               , the standard deviation of the Gaussian
%                                filter in pixels
% * 'guessssfactor'            , the integer specifying the degree of data
%                                subsampling
% * 'guessnBins'               , number of radial bins to determine the
%                                point cloud in
% * 'guessrmRadialOutliers'    , the threshold to remove radial outliers
% * 'guessrmIntensityOutliers' , the threshold to remove intensity outliers
%

guessChannel             = 1 ;
guessSigma               = 2 ;
guessssfactor            = 8 ;
guessnBins               = 120 ;
guessrmRadialOutliers    = 1.2 ;
guessrmIntensityOutliers = 1.2 ;

%
% Options for visualizing the intial point cloud.
%
% * 'plot2Ddimval'   , the dimension along which to take cross section
% * 'plot2Dval'      , the slice to use for plotting
% * 'plot3Dssfactor' , the factor to reduce the number of points shown
%

plot2Ddimval   = 'z' ;
plot2Dval      = 350 ;
plot3Dssfactor = 4 ;

%
% Options for coarse fitting the surface.
%
% * 'Rfit'     , the degree of polynomial fit for radius
% * 'Xfit'     , the degree of polynomial fit for x
% * 'Yfit'     , the degree of polynomial fit for y
% * 'efit'     , the degree of polynomial fit for ellipse eccentricity
% * 'phasefit' , the degree of polynomial fit for ellipse orientation
%

Rfit     = 6 ;
Xfit     = 4 ;
Yfit     = 4 ;
efit     = 3 ;
phasefit = 0 ;

%
% Options for morphsnakes surface detection.
%
% * 'ms_scriptDir'      , the location of the morphsnakes repo
% * 'channel'           , 
% * 'foreGroundChannel' , 
% * 'niter'             , 
% * 'niter0'            , 
% * 'lambda1'           , 
% * 'lambda2'           , 
% * 'nu'                , 
% * 'smoothing'         , 
% * 'post_nu'           , 
% * 'post_smoothing'    , 
% * 'exit_thres'        , 
% * 'pre_nu'            ,
% * 'pre_smoothing'     , 
% * 'radius_guess'      , 
% * 'clip'              , the boolean wheter to clip the intensity of the
%                         raw data
% * 'save'              , the boolean whether to save intermediate steps
% * 'mlxprogram'        , the location of the mesh script to use
%

ms_scriptDir      = '/mnt/data/code/morphsnakes_wrapper/' ;
channel           = -1 ;
foreGroundChannel = 1 ;
niter             = 40 ;
niter0            = 40 ;
lambda1           = 1 ;
lambda2           = 2 ;
nu                = 0.1 ;
smoothing         = 1 ;
post_nu           = -1 ;
post_smoothing    = 1 ;
exit_thres        = 0.0001 ;
pre_nu            = 0 ;
pre_smoothing     = 0 ;
radius_guess      = 30 ;
clip              = 800.0 ;
save              = false ;
mlxprogram        = '/mnt/data/code/meshlab_codes/surface_rm_resample2k_reconstruct_LS3_ssfactor4.mlx' ;

%
% Options for loading the mesh from morphsnakes.
%
% * 'normal_step_for_mesh' , the number of pixels to nomally evolve
%

normal_step_for_mesh = 0 ;

%
% Options for setting the orientation of the system.
%
% * 'ventralpresent' , the boolean whether the ventral furrow is present
% * 'celllayer'      , the underestimated guess of the thickness of the
%                      embryo cell layer in pixels
% * 'celllayerstep'  , the number of pixels to increase cell layer guess by
%                      each iteration
% * 'min_thres_wo'   , the minimum threshold for x and y eigenvalues to use
%                      if the ventral furrow is not present
%

ventralpresent = false ;
celllayer      = 30 ;
celllayerstep  = 7 ;
min_thres_wo   = 0.1 ;

%
% Options for inspecting morphsnakes fit.
%
% * 'mesh2Ddim'     , the dimension along which to visualize
% * 'mesh2Dzval'    , the value of the slice to visualize
% * 'mesh2Dnoalign' , 
%

mesh2Ddim     = 'z' ;
mesh2Dzval    = 350 ;
mesh2Dnoalign = true ;

%
% Options for the pullback styles.
%
% * 'cylinderone'     , the boolean whether to compute cylinder one
% * 'cylindertwo'     , the boolean whether to compute cylinder two
% * 'cylinderoneprop' , the boolean whether to compute cylinder one proper
% * 'cylindertwoprop' , the boolean whether to compute cylinder two proper

cylinderone     = true ;
cylindertwo     = true ;
cylinderoneprop = true ;
cylindertwoprop = true ;

%
% Options for the onion layers.
%
% * 'onionnLayers'       , the number of onion layers is an odd integer
% * 'onionlayerDistance' , the onion layer distance in pixels
% * 'onionsigma'         , the sigma used for smoothing onion layers
% * 'onionmakeIP'        , the options are SIP, MIP, or both
% * 'onionIPonly'        , 
%

onionnLayers       = 3 ;
onionlayerDistance = 5 ;
onionsigma         = 10 ;
onionmakeIP        = 'both' ;
onionIPonly        = false ;

%
% Options for saving to disc.
%
% * 'soiDir'         , the directory to save the pullbacks to
% * 'imwriteOptions' , the image writing options
% * 'make8bit'       , the boolean whether to save files as 8bit
%

soiDir         = 'cylinder_soi_t%d_c%d' ;
imwriteOptions = {'tif'} ;
make8bit       = false ;

%% Automatic Pullback Generation
%
% Automatically generate pullbacks for all times and channels.
%

if detectchan
    % Detect the number of channels in the folders
    filenamewild = [fnfull '.tif'] ;
    filenamewild = sprintf(filenamewild, firstTime) ;
    filenamewild = strrep(filenamewild, 'c', 'c*') ;
    numChannels = length(dir(filenamewild)) ;
else
    numChannels = UserChan ;
end

if detecttime
    % Detect the number of time points in each channel
    filenamewild = strrep([fnfull '.tif'], 'Ch', ['Ch' num2str(firstChan)]) ;
    filenamewild = strrep(filenamewild, 'TP%d', 'TP*') ;
    numTimes = length(dir(filenamewild)) ;
else
    numTimes = UserTimes ;
end

% Set the number of time points available
timePointsava   = firstTime : (firstTime + numTimes - 1) ;

for c = firstChan : (firstChan + numChannels - 1)
    
    fn = strrep(fnfull, 'c', ['c' num2str(c)]) ;
    
    for t = firstTime : (firstTime + numTimes - 1)
        
        projectDir = cd ;
        dataDir = cd ;
        
        parametersfile = ['pipeparameters_c' num2str(c) '_t' num2str(t) '.txt'] ;

        if ~isfile(parametersfile) || overwrite

            %
            % Start by creating an experiment object, optionally pass on
            % the project directory (otherwise it will ask). This serves as
            % a frontend for data loading, detection, fitting, etc.
            %

            xp = project.Experiment(projectDir, dataDir) ;

            %
            % Set required additional information on the files.
            %
            % We assume one individual image stack for each time point,
            % labeled by time. To be able to load the stack, we need to
            % tell the project where the data is, what convention is
            % assumed for the file names, available time points and the
            % stack resolution. Options to modules in ImSAnE are organised
            % in matlab structures, that is a pair of field name and value
            % are provided for each option.
            %

            fileMeta                 = struct() ;
            fileMeta.dataDir         = dataDir ;
            fileMeta.filenameFormat  = [fn '.tif'] ;
            fileMeta.timePoints      = timePointsava ;
            fileMeta.stackResolution = stackResolution ;
            fileMeta.swapZT          = swapZT ;

            %
            % Set required additional information on the experiment. A verbal data set
            % description, Jitter correct by translating the sample, which time point to 
            % use for fitting, etc.
            %
            % The following project metadata information is required 
            %
            % * 'channelsUsed'   , the channels used, e.g. [1 3] for RGB
            % * 'channelColor'   , mapping from element in channels used to RGB = 123
            %

            expMeta                  = struct() ;
            expMeta.channelsUsed     = channelsUsed ;
            expMeta.channelColor     = channelColor ;
            expMeta.description      = description ;
            expMeta.dynamicSurface   = 0 ;
            expMeta.jitterCorrection = 0 ;
            expMeta.fitTime          = fileMeta.timePoints(t - firstTime + 1) ; 
            expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector' ;
            expMeta.fitterType       = 'surfaceFitting.spherelikeFitter' ;

            %
            % Now set the meta data in the experiment.
            %

            xp.setFileMeta(fileMeta) ;
            xp.setExpMeta(expMeta) ;

            %
            % initNew() reads the stack size from the first available time
            % point, then initializes the fitter and deterctor and creates
            % fitOptions and detectOptions based on their defaults.
            %

            xp.initNew() ;

            %
            % loadTime sets currentTime, loads the stack, and resets the
            % detector and fitter with the appropriate options for that
            % time. Optionally rescale the stack to unit aspect ratio.
            %

            xp.loadTime(t) ;
            xp.rescaleStackToUnitAspect() ;

            %
            % Create the struct of detect options.
            %

            myDetectOpts = struct('channel', guessChannel, ...
                'sigma', guessSigma, ...
                'ssfactor', guessssfactor, ...
                'nBins', guessnBins, ...
                'rmRadialOutliers', guessrmRadialOutliers, ...
                'rmIntensityOutliers', guessrmIntensityOutliers, ...
                'zDim', zDim) ;

            %
            % Set the detect options in the project
            %

            xp.setDetectOptions(myDetectOpts) ;

            %
            % Calling detectSurface runs the surface detector and creates
            % the detector.pointCloud object.
            %

            xp.detectSurface() ;

            %
            % Inspect the point cloud over a cross section in the data.
            %

            if visualize
                inspectOptions = struct('dimension', plot2Ddimval, ...
                    'value', plot2Dval, ...
                    'pointCloud', 'r') ;
                figure ;
                xp.detector.inspectQuality(inspectOptions, xp.stack) ;
            end


            %
            % Plot the points on the surface as a 3D point cloud.
            %

            if visualize
                xp.detector.pointCloud.inspect(plot3Dssfactor) ;
            end

            %
            % Fit the surface coarsly to prepare levelset.
            %
            % We fit the surface in a cylindrical basis, with radius,
            % eccentricity, centre of mass, and ellipse orientation as
            % slowly varying polynomials of the z axis.
            %
            % The fitter behaves similar to the detector, but the fitting
            % of the surface is done by calling the function fitSurface()
            % in the experiment class. This function calls
            % fitter.fitSurface(detector.pointCloud), but also stores the
            % resulting fittedParam in the experiment class.
            %
            % Initially we fit the surface with a coarse set of values, to
            % obtain the initial levelset for morphsnakes.
            %

            fitOptions = struct('R', Rfit, ...
                'X', Xfit, ...
                'Y', Yfit, ...
                'e', efit, ...
                'phase', phasefit, ...
                'path', fullfile(projectDir, 'debugOutput')) ;
            xp.setFitOptions(fitOptions) ;
            xp.fitSurface() ;

            %
            % Now use the polynomial fit of z of ellipses to get inside/outside
            %
            % If initial guess does not exist, create it.
            %

            initls_h5fn = sprintf('Time_%06d_c%d_levelset.h5', xp.currentTime, c) ;
            if ~exist(initls_h5fn, 'file')
                % Create initial level set from fit
                X0 = @(Z) polyval(xp.fittedParam.pX, double(Z), xp.fittedParam.SX, xp.fittedParam.muX) ;
                Y0 = @(Z) polyval(xp.fittedParam.pY, double(Z), xp.fittedParam.SY, xp.fittedParam.muY) ;
                rsq = @(Z) polyval(xp.fittedParam.pR, double(Z), xp.fittedParam.SR, xp.fittedParam.muR) ;
                ee = @(Z) polyval(xp.fittedParam.pe, double(Z), xp.fittedParam.Se, xp.fittedParam.mue) ;
                % phi = @(Z) polyval(xp.fittedParam.pphase, double(Z), xp.fittedParam.Sphase, xp.fittedParam.muphase) ;
                % Build the level set from the subsampled data.
                opts = xp.detector.options ; 
                data = xp.stack.image.apply{opts.channel} ;
                data = data(1:opts.ssfactor:end, 1:opts.ssfactor:end, 1:opts.ssfactor:end) ;

                ls = 0 * data ;
                zvals = (1:size(data, zDim)) ;
                if zDim == 1
                    xyinds = [3, 2] ;
                elseif zDim == 2
                    error('Decide the order of the axes here.')
                end
                xarr = (1:size(data, xyinds(1))) * opts.ssfactor ;
                yarr = (1:size(data, xyinds(2))) * opts.ssfactor ;
                [xgrid, ygrid] = meshgrid(xarr, yarr) ;
                for zind = zvals
                    zval = zind * opts.ssfactor ;
                    if mod(zval, 100) == 0
                        disp(['zval = ', num2str(zval)])
                    end
                    cosp = cos(0) ;
                    sinp = sin(0) ;
                    % If phase of the ellipse is necessary uncomment below.
                    %
                    % cosp = cos(phi(zval)) ;
                    % sinp = sin(phi(zval)) ;
                    %
                    % Get the xy values inside the polynomial.
                    rsqz = rsq(zval) ;
                    rbigsqz = rsqz / (1 - ee(zval)^2 )^2 ;
                    X0z = X0(zval) ;
                    Y0z = Y0(zval) ; 
                    dd = (((xgrid - X0z) * cosp + (ygrid - Y0z) * sinp).^2 / rsqz + ((xgrid - X0z) * sinp - (ygrid - Y0z) * cosp).^2 / rbigsqz) ;
                    ls(zind,:,:) = (dd > 0) & (dd < 1) ;
                end
                % Write the levelset to disk.
                h5create(initls_h5fn,'/implicit_levelset', size(ls))
                h5write(initls_h5fn, '/implicit_levelset', ls)
            else
                disp('implicit_levelset found on disk')
            end

            %
            % Initialize morphsnakes detectors.
            %

            expMeta.detectorType = 'surfaceDetection.integralDetector_rawdata' ;
            expMeta.fitterType = 'surfaceFitting.cylinderMeshWrapper' ;
            xp.setExpMeta(expMeta) ;
            xp.initNew()

            %
            % If you initialize a new project with xp.initNew(),
            % xp.currentTime returns 1 rather than the time point that we
            % are working with, which is why we ned to use xp.setTime().
            %

            xp.setTime(t) ;

            %
            % File naming for morphsnakes.
            %

            ssfactor = xp.detector.options.ssfactor ;
            ofn_ply = ['mesh_apical_ms_t' num2str(t) '_c' num2str(c)] ;
            ofn_ls = ['msls_apical_t' num2str(t) '_c' num2str(c)] ;
            ofn_smoothply = ['mesh_apical_t' num2str(t) '_c' num2str(c)] ;

            %
            % Name the output mesh directory
            %

            msls_exten = ['_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ] ;
            msls_exten = [msls_exten '_s' num2str(smoothing, '%d') ] ;
            msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps', num2str(post_smoothing)] ;
            msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ] ;
            if projectDir(end) ~= '/'
                projectDir = [projectDir '/'] ;
            end
            mslsDir = [projectDir 'msls_output_t' num2str(t) '_c' num2str(c)] ;
            mslsDir = [mslsDir msls_exten '/'] ;
            if ~exist(mslsDir,'dir')
                mkdir(mslsDir)
            end

            %
            % Set the detection options for morphsnakes.
            %

            msDetectOpts = struct('channel', channel, ...
                'ssfactor', ssfactor, ...
                'foreGroundChannel', foreGroundChannel, ...
                'niter', niter, ...
                'niter0', niter0, ...
                'lambda1', lambda1, ...
                'lambda2', lambda2, ...
                'nu', nu, ...
                'smoothing', smoothing, ...
                'post_nu', post_nu, ...
                'post_smoothing', post_smoothing, ...
                'exit_thres', exit_thres, ...
                'fileName', sprintf(fn, xp.currentTime), ...
                'mslsDir', mslsDir, ...
                'ofn_ls', ofn_ls, ...
                'ofn_ply', ofn_ply,...
                'ms_scriptDir', ms_scriptDir, ...
                'timepoint', xp.currentTime, ...
                'zdim', zDim, ...
                'pre_nu', pre_nu, ...
                'pre_smoothing', pre_smoothing, ...
                'ofn_smoothply', ofn_smoothply, ...
                'mlxprogram', mlxprogram, ...
                'init_ls_fn', initls_h5fn, ...
                'run_full_dataset', 'none', ...
                'radius_guess', radius_guess, ...
                'dset_name', 'inputData', ...
                'clip', clip, ...
                'save', save, ...
                'center_guess', 'empty_string', ...
                'plot_mesh3d', false, ...
                'mask', 'none') ;

            %
            % Set the morphsnakes detect options.
            %

            xp.setDetectOptions(msDetectOpts) ;

            %
            % Move the initial level set to the msls dir.
            %

            if exist(initls_h5fn, 'file')
                system(['mv ' initls_h5fn ' ' fullfile(mslsDir, initls_h5fn)]) ;
            end

            %
            % Create downsampled dataset as hdf5
            %

            if ~exist([projectDir sprintf(fn, xp.currentTime) '.h5'], 'file')
                xp.detector.prepareDownsampled(xp.stack) ;
                disp('Done outputting downsampled data h5 for surface detection.')
            else
                disp('The h5 file was already output.')
            end

            disp('Downsampling to h5 is done.')

            %
            % Detect the surface of the embryo.
            %

            xp.detectSurface()

            %
            % Read in the mesh file.
            %

            mesh_outfn = [ofn_smoothply, num2str(fileMeta.timePoints(xp.currentTime - firstTime + 1), '%06d'), '.ply'] ;
            outputMesh = fullfile(mslsDir, mesh_outfn) ;
            mesh = read_ply_mod(outputMesh) ;
            msls_axis_order = 'zxyc' ;

            %
            % There is a rotation implicit here by permuting axes
            %

            if strcmp(msls_axis_order, 'zxyc')
                mesh.v  = mesh.v(:, [2,3,1]) * ssfactor ;
                mesh.vn = mesh.vn(:, [2,3,1]) * ssfactor ;
            end

            %
            % Make sure vertex normals are normalized.
            %

            mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) ) ;

            %
            % Normally evolve vertices.
            %

            mesh.v = mesh.v + normal_step_for_mesh .* mesh.vn ;
            
            if ventralpresent
                %
                % File naming conventions.
                %

                datafilename = [sprintf(fn, xp.currentTime) '.h5'] ;
                data = h5read(datafilename, '/inputData') ;
                levelset = h5read(fullfile(mslsDir, initls_h5fn), '/implicit_levelset') ;

                %
                % Determine what slice to use.
                %

                S = sum(levelset, 2);
                S = squeeze(S) ;
                S = sum(S, 2) ;
                S = squeeze(S) ;
                slice = round(mean([find(S,1) find(S,1,'last')]));

                %
                % Create binary image from particular slice of data.
                %

                imageslice = data(slice, :, :) ;
                imageslice = squeeze(imageslice) ;
                imageslice = kron(imageslice, ones(opts.ssfactor)) ;
                imageslicebw = imbinarize(imageslice, 550) ;

                CC = bwconncomp(imageslicebw) ;

                %
                % Clean up the image of extra regions that happen to be above
                % the threshold.
                %

                numobjects = CC.NumObjects ;
                if numobjects > 1
                    for a = 1 : numobjects - 1
                        numPixels = cellfun(@numel,CC.PixelIdxList) ;
                        [biggest,idx] = min(numPixels) ;
                        imageslicebw(CC.PixelIdxList{idx}) = 0 ;
                        CC = bwconncomp(imageslicebw) ;
                    end
                end

                %
                % Keep an extra copy of clean image.
                %

                imageslicebwnofill = imageslicebw ;

                %
                % Fill image for center of region calculation.
                %

                imageslicebwfill = imfill(imageslicebw,'holes') ;

                %
                % Erode non-filled image of the embryo to leave the ventral
                % furrow and clean up the image.
                %

                celllayer = 30 ;
                se = strel('disk', celllayer) ;
                imageslicebw = imerode(imageslicebwnofill, se) ;
                CC = bwconncomp(imageslicebw) ;
                numobjects = CC.NumObjects ;

                while numobjects < 2
                    se = strel('disk', celllayer) ;
                    imageslicebw = imerode(imageslicebwnofill, se) ;
                    CC = bwconncomp(imageslicebw) ;
                    numobjects = CC.NumObjects ;
                    celllayer = celllayer + 5 ;
                end

                if numobjects > 1
                    for a = 1 : numobjects - 1
                        numPixels = cellfun(@numel,CC.PixelIdxList) ;
                        [biggest,idx] = min(numPixels) ;
                        imageslicebw(CC.PixelIdxList{idx}) = 0 ;
                        CC = bwconncomp(imageslicebw) ;
                    end
                end

                center = regionprops(imageslicebwfill, 'Centroid') ;
                centernofill = regionprops(imageslicebw, 'Centroid') ;

                %
                % Plot the line from the center of mass of the embryo to the
                % center of mass of the ventral furrow.
                %

                if visualize
                    figure ;
                    imagesc(imageslicebwnofill)
                    hold on ;
                    plot(center.Centroid(1), center.Centroid(2), 'go')
                    plot(centernofill.Centroid(1), centernofill.Centroid(2), 'ro')
                    plot([center.Centroid(1) centernofill.Centroid(1)], [center.Centroid(2) centernofill.Centroid(2)], 'r')
                end

                %
                % Compute the phase of the DV axis
                %

                phase = atan2((centernofill.Centroid(2) - center.Centroid(2)), (centernofill.Centroid(1) - center.Centroid(1))) ;
                min_thres = 0.1 ;
                
            else
                
                phase = 0 ;
                min_thres = min_thres_wo ;
                
            end
            
            %
            % Determine the orientation of the region of interest.
            %
            % Note that the first argument of pc.determineROI is the
            % margin around the point cloud.
            %
            points    = mesh.v;
            pc        = surfaceDetection.PointCloud(points);
            pc.determineROI(1, min_thres) ;
            translation = pc.ROI.translation ;
            rotation    = pc.ROI.rotation ;

            %
            % Set the fit options for the surfacefitter.
            %
            fitOptions = struct('chartSeeds', [], ...
                'diskSeeds', [], ...
                'phase', phase, ...
                'transitionWidth', 0, ...
                'fixAxis', 1, ...
                'rotation', rotation, ...
                'translation', translation, ...
                'fixResolution', 1, ...
                'resolution', []) ;
            xp.setFitOptions(fitOptions) ;
            xp.fitSurface(mesh) ;

            %
            % Inspect fit and point cloud over a cross section in the data.
            %
            if visualize
                inspectOptions = struct('dimension', mesh2Ddim, ...
                    'value', mesh2Dzval, ...
                    'pointCloud', 'b', ...
                    'noalign', mesh2Dnoalign) ;
                figure ;
                xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack) ;
            end

            %
            % Inspect the whole mesh in three-dimensions. Note that the y
            % axis should be the AP axis of the embryo.
            %
            if visualize
                figure ;
                xp.fitter.inspectMesh() ;
                xlabel('x') ;
                ylabel('y') ;
                zlabel('z') ;
            end

            %
            % Define the pullback styles.
            %
            xp.fitter.setDesiredChart('cylinder1', cylinderone) ;
            xp.fitter.setDesiredChart('cylinder2', cylindertwo) ;
            xp.fitter.setDesiredChart('cylinder1_proper', cylinderoneprop) ;
            xp.fitter.setDesiredChart('cylinder2_proper', cylindertwoprop) ;
            xp.generateSOI() ;
            disp('Done generating SOI for first timepoint.')

            %
            % Set the options for the onion projections.
            %
            onionOpts = struct('nLayers', onionnLayers, ...
                'layerDistance', onionlayerDistance, ...
                'sigma', onionsigma, ...
                'makeIP', onionmakeIP, ...
                'IPonly', onionIPonly) ;

            %
            % Pass the region of interest and the current time to pull back
            % the stack in the desired charts. This generates the data
            % fields containing the pullback.
            %
            xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts) ;
            disp('Done pulling back onions for first timepoint.')

            %
            % Now we extract the data field from the surface of interest at
            % the current time, which is the time of the fit.
            %
            fitOptions    = xp.fitter.fitOptions ;
            data          = xp.SOI.getField('data_MIP') ;
            data          = data(xp.tIdx(xp.currentTime)) ;
            type          = 'cylinder' ;
            patchName     = 'cylinder2_index' ;
            transformName = 'cylinder2' ;
            pb = data.getPatch(patchName).getTransform(transformName).apply{1} ;
            if visualize
                figure, imshow(pb',[],'InitialMagnification',66) ;
                figure ; imagesc(pb') ; colorbar ;
            end

            %
            % Save the pullbacks to disc.
            %
            saveDir = fullfile(projectDir, sprintf(soiDir, xp.currentTime, c)) ;
            options = struct('dir', saveDir, ...
                'imwriteOptions', {imwriteOptions}, ...
                'make8bit', make8bit) ;
            xp.SOI.save(options)

            %
            % Save a plain text file containing the options used in the
            % experiemt.
            %
            optionstring = ['guessChannel = ' num2str(guessChannel), ...
                '\nguessSigma = ' num2str(guessSigma), ...
                '\nguessssfactor = ' num2str(guessssfactor), ...
                '\nguessnBins = ' num2str(guessnBins), ...
                '\nguessrmRadialOutliers = ' num2str(guessrmRadialOutliers), ...
                '\nguessrmIntensityOutliers = ' num2str(guessrmIntensityOutliers), ...
                '\nzDim = ' num2str(zDim), ...
                '\nRfit = ' num2str(Rfit), ...
                '\nXfit = ' num2str(Xfit), ...
                '\nYfit = ' num2str(Yfit), ...
                '\nefit = ' num2str(efit), ...
                '\nphasefit = ' num2str(phasefit), ...
                '\nchannel = ' num2str(channel), ...
                '\nforeGroundChannel = ' num2str(foreGroundChannel), ...
                '\nniter = ' num2str(niter), ...
                '\nniter0 = ' num2str(niter0), ...
                '\nmlx_program = ' mlxprogram, ...
                '\nlambda1 = ' num2str(lambda1), ...
                '\nlambda2 = ' num2str(lambda2), ...
                '\nexit_thres = ' num2str(exit_thres), ...
                '\nsmoothing = ' num2str(smoothing), ...
                '\nnu = ' num2str(nu), ...
                '\npre_nu = ' num2str(pre_nu), ...
                '\npre_smoothing = ' num2str(pre_smoothing), ...
                '\npost_nu = ' num2str(post_nu), ...
                '\npost_smoothing = ' num2str(post_smoothing), ...
                '\nradius_guess = ' num2str(radius_guess), ...
                '\nclip = ' num2str(clip), ...
                '\nonionnLayers = ' num2str(onionnLayers), ...
                '\nonionlayerDistance = ' num2str(onionlayerDistance), ...
                '\nonionsigma = ' num2str(onionsigma)] ;

            fparam = fopen(parametersfile, 'wt') ;
            fprintf(fparam, optionstring) ;
            fclose(fparam) ;
            
        else
            
            disp(['Data in ' projectDir ' has already been processed.'])
            
        end
     
    end
    
end
