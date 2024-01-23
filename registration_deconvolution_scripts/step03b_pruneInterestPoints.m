% Remove interest points based on local density (determined via kdTree)
% NPMitchell 2023


rootdir = './interestpoints/' ;
tps = 0:101 ;   % timepoints to prune
vtiles = 0:15 ; % view tiles whose interest points to prune

% Parameters for the experimental data:
dz = 1.4 ; % um, step size between frames in z
dx = 0.2619 ; % um, step size between pixels in xy

% Parameters for the pruning based on Density
useDensity = true ;
densityThres = 0.02 ; % threshold local density
preview = false ; % show pruning in figure window
pausetime = 0.1 ; % time to show pruning, in seconds

% Parameters for the pruning based on ROI
% Store limits and radius (if cylinder) of ROI to include/exclude
useROI = false ;
ROIs = {[]} ; % if useROI is true, define ROIs here as cell of arrays
ROIstyle = 'cylinder' ; % 'box' (xmin, xmax, ymin, ymax, zmin, zmax)
                        % or 'cylinder' (x,y,radius, zmin, zmax)

% include or exclude interestpoints in the ROI / above density threshold
includeExclude = 'exclude' ; 

% Preallocate
clc ; 
colors = define_colors(length(vtiles)) ;
% Count number removed
nipsRm = repmat({zeros(length(tps), 1)},length(vtiles), 1); 
nipsKeep = nipsRm ;

for tidx = 1:length(tps)
    tp = tps(tidx) ;
    if mod(tp, 10) == 0
        disp(['pruning interest points for tp=' num2str(tp) '/' num2str(max(tps))])
    end
    for vId = 1:length(vtiles) 
        fn2find = fullfile(rootdir, sprintf('tpId_%d_viewSetupId_%d.beads.ip.txt', tp, vtiles(vId))) ;
        fns = dir(fn2find) ;
        if length(fns) == 1
            fn = fullfile(fns(1).folder, fns(1).name) ;
            % disp([' > inspecting ' fn])
            ips = dlmread(fn, '\t', 1, 0) ;
            assert(size(ips, 1) == max(ips(:, 1))+1) ; % check that the #ips is correct
            ips = ips(:, 2:end) ;

            % Prune IPs
            % scatter3(ips(:, 1), ips(:, 2), ips(:, 3), 10, 'filled')
            % xlabel('x'); ylabel('y'), zlabel('z'); axis equal
            
            if useDensity
                % rescale pointCloud based on dz
                pointCloud = ips ;
                pointCloud(:, 3) = pointCloud(:, 3) * dz/dx ;
    
                % Create a KD-Tree for fast nearest neighbor search
                kdtree = KDTreeSearcher(pointCloud);
                
                % Use knnsearch to find the indices of the k-nearest neighbors and their distances
                [indices, distances] = knnsearch(kdtree, pointCloud, 'K', 10);
    
                % Calculate the local density based on the distances to the k-th nearest neighbor
                localDensities = 1 ./ mean(distances, 2);% Scale local densities to fit the colormap
                cmin = min(localDensities);
                cmax = max(localDensities);
                scaledDensities = (localDensities - cmin) / (cmax - cmin);
                
                % Now, the scaledDensities array contains values between 0 and 1 based on local densities.
                if strcmpi(includeExclude, 'exclude')
                    reject = localDensities > densityThres ;
                elseif strcmpi(includeExclude, 'include')
                    reject = localDensities < densityThres ;
                else
                    error('Must set includeExclude to <include> or <exclude>')
                end
                % Plot the point cloud with colors determined by local densities
                if preview
                    clf
                    colormap bwr
                    subplot(1, 2, 1)
                    scatter3(ips(:, 1), ips(:, 2), ips(:, 3), 10, localDensities, 'filled');
                    xlabel('x'); ylabel('y'), zlabel('z'); axis equal
                    caxis([0, densityThres])
                    view(2) ; colorbar('south')
                    subplot(1, 2, 2)
                    scatter3(ips(:, 1), ips(:, 2), ips(:, 3), 10, reject, 'filled');
                    xlabel('x'); ylabel('y'), zlabel('z'); axis equal
                    view(2)
                    sgtitle(['tp=' num2str(tp) ': view ' num2str(vtiles(vId))])
                    set(gcf, 'color', 'w')
                    pause(pausetime)
                end
            end
            
            if useROI
                error('code up the ROI step here, independent of whether we already pruned based on density')
            end
            
            % Create output interest point matrix
            ips = ips(~reject, :) ;
            newIds = 0:(size(ips, 1)-1) ;
            outputIP = [newIds', ips(:, 1), ips(:, 2), ips(:, 3)] ;

            % Custom write txt to disk
            fid = fopen(fn, 'wt');
            % make header
            header = 'id	x	y	z' ;
            fprintf(fid, [header '\n']);  
            for ii = newIds
                fprintf(fid, [num2str(ii), '\t', ...
                    sprintf('%0.14f', ips(ii+1, 1)), '\t', ...
                    sprintf('%0.14f', ips(ii+1, 2)), '\t', ...
                    sprintf('%0.14f', ips(ii+1, 3)), '\n' ]) ;
            end
            fclose(fid);

            % This doesn't work
            % write_txt_with_header(fn, outputIP, header, '\t', '%0.14f')
            nipsRm{vId}(tidx) = sum(reject) ;
            nipsKeep{vId}(tidx) = size(ips, 1) ;
        else
            msg = 'Found more than one interestpoint file for this timepoint/viewtile: tp=%d vtile=%d';
            error(sprintf(msg, tp, vtile))
        end
    end
end


%% Plot the # interest points as a function of time for each view
close all; 
legEntry = cell(1, length(vtiles)) ;
for vId = 1:length(vtiles)
    plot(tps, nipsKeep{vId}, '.-', 'color', colors(vId, :))
    legEntry{vId} = sprintf('View %d', vtiles(vId)) ;
    hold on;
end
xlabel('timepoint')
ylabel('#interest points')
legend(legEntry, 'location', 'eastOutside')
title('Counts after interest point pruning')
currentDateTime = datetime('now', 'Format', 'yyyyMMddHHmm');
dateTimeString = char(currentDateTime);

saveas(gcf, ['interestpoint_counts_Pruned_' dateTimeString '.png'])
