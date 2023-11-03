% Evaluate # of interestpoints per view

rootdir = './interestpoints/' ;
tps = 0:101 ;   % timepoints
vtiles = 8:15 ; % view tiles
previewIPs = false ; 
pauseTime = 0.1 ;

% Preallocate
clc ; 
colors = define_colors(length(vtiles)) ;
nips = repmat({zeros(length(tps), 1)},length(vtiles), 1); 

for tidx = 1:length(tps)
    tp = tps(tidx) ;
    if mod(tp, 10) == 0
        disp(['counting interest points for tp=' num2str(tp) '/' num2str(max(tps))])
    end
    for vId = 1:length(vtiles) 
        fn2find = fullfile(rootdir, sprintf('tpId_%d_viewSetupId_%d.beads.ip.txt', tp, vtiles(vId))) ;
        fns = dir(fn2find) ;
        if length(fns) == 1
            fn = fullfile(fns(1).folder, fns(1).name) ;
            % disp([' > inspecting ' fn])
            ips = dlmread(fn, '\t', 1, 0) ;
            nips{vId}(tidx) = size(ips, 1) ;
            assert(size(ips, 1) == max(ips(:, 1))+1) ; % check that the #ips is correct
        
        
            if previewIPs
                % rescale pointCloud based on dz
                ips = ips(:, 2:end) ;
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
    
                % Plot the point cloud with colors determined by local densities
                clf
                colormap bwr
                scatter3(ips(:, 1), ips(:, 2), ips(:, 3), 10, localDensities, 'filled');
                xlabel('x'); ylabel('y'), zlabel('z'); axis equal
                colorbar
                title(['tp=' num2str(tp) ': view ' num2str(vtiles(vId))])
                pause(pauseTime)
            end
        
        
        else
            msg = 'Found more/less than one interestpoint file for this timepoint/viewtile: tp=%d vtile=%d';
            error(sprintf(msg, tp, vtile))
        end
    end
end

%% Plot the # interest points as a function of time for each view
close all; 
legEntry = cell(1, length(vtiles)) ;
for vId = 1:length(vtiles)
    plot(tps, nips{vId}, '.-', 'color', colors(vId, :))
    legEntry{vId} = sprintf('View %d', vtiles(vId)) ;
    hold on;
end
xlabel('timepoint')
ylabel('#interest points')
legend(legEntry)
currentDateTime = datetime('now', 'Format', 'yyyyMMddHHmm');
dateTimeString = char(currentDateTime);
saveas(gcf, ['interestpoint_counts_' dateTimeString '.png'])


