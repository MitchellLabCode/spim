% Evaluate # of interestpoints per view
close all 
clear all

rootdir = './interestpoints/' ;
previewIPs = true ; 
pauseTime = 0.1 ;
previewEveryN = 10 ;
densityThres = 0.02 ; % threshold local density
dz = 1.0 ; % um / page, axial sampling
dx = 0.195 ; % um / pixel, in-plane sampling

%% Automatically determine timepoints and view tiles
% Automatically determine timepoints and view tiles
tp_files = dir('./c*_t*_a*.ome.tif');
tps = unique(cellfun(@(x) sscanf(x, 'c%*d_t%d_a%*d.ome.tif'), {tp_files.name}));

% Extract the channel and angle info
vtiles = cellfun(@(x) sscanf(x, 'c%d_t%*d_a%d.ome.tif', [1, 2]), {tp_files.name}, 'UniformOutput', false);

% Convert vtiles into a matrix for unique operation
vtiles_matrix = cell2mat(vtiles');  % Convert to 2D numeric matrix (channel, angle)
vtiles = unique(vtiles_matrix, 'rows');  % Get unique rows (tuples of channel, angle)
vtiles = 0:size(vtiles, 1)-1 ;

%% Count interestpoints
% Preallocate
clc ; 
colors = lines(length(vtiles)) ;
nips = repmat({zeros(length(tps), 1)},length(vtiles), 1); 
missing_files = {};

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
            ips = dlmread(fn, '\t', 1, 0) ;
            nips{vId}(tidx) = size(ips, 1) ;
            assert(size(ips, 1) == max(ips(:, 1))+1) ;

            if previewIPs && mod(tidx, previewEveryN) == 1
                ips = ips(:, 2:end) ;
                pointCloud = ips ;
                pointCloud(:, 3) = pointCloud(:, 3) * dz/dx ;
                
                kdtree = KDTreeSearcher(pointCloud);
                [indices, distances] = knnsearch(kdtree, pointCloud, 'K', 10);
                localDensities = 1 ./ mean(distances, 2);
                cmin = min(localDensities);
                cmax = max(localDensities);
                scaledDensities = (localDensities - cmin) / (cmax - cmin);

                clf
                colormap bwr
                scatter3(ips(:, 1), ips(:, 2), ips(:, 3), 10, localDensities, 'filled');
                xlabel('x'); ylabel('y'), zlabel('z'); axis equal
                caxis([0, densityThres])
                view(2)
                colorbar
                title(['tp=' num2str(tp) ': view ' num2str(vtiles(vId))])
                pause(pauseTime)
            end
        else
            missing_files{end+1} = sprintf('tp=%d vtile=%d', tp, vtiles(vId));
        end
    end
end

% Display missing files
if ~isempty(missing_files)
    disp('Missing files:')
    fprintf('%s\n', missing_files{:});
else
    disp('All expected files are present.')
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
legend(legEntry, 'location', 'northeastoutside')
currentDateTime = datetime('now', 'Format', 'yyyyMMddHHmm');
dateTimeString = char(currentDateTime);
set(gcf, 'color', 'w')
title('Interest points summary')
saveas(gcf, ['interestpoint_counts_' dateTimeString '.png'])
