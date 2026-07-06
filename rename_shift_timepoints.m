function rename_shift_timepoints(masterDir, desiredFirstTimepoint)
%SHIFT_TIMEPOINTS Shift Cam_*_#####.* timepoints inside stack_* folders.
%
% Usage:
%   shift_timepoints('/path/to/part2/raw', 21)
%
% This makes old 00000 -> 00021, old 00001 -> 00022, etc.
%
% It searches recursively for folders starting with "stack_" and renames
% files matching: Cam_<something>_00000.<ext> (e.g. .json, .lux.h5).

    if nargin < 2
        error('Usage: shift_timepoints(masterDir, desiredFirstTimepoint)');
    end
    if ~isfolder(masterDir)
        error('masterDir is not a folder: %s', masterDir);
    end
    if ~isscalar(desiredFirstTimepoint) || desiredFirstTimepoint < 0 || fix(desiredFirstTimepoint) ~= desiredFirstTimepoint
        error('desiredFirstTimepoint must be a nonnegative integer.');
    end

    % Find all directories under masterDir
    allDirs = dir(fullfile(masterDir, '**'));
    allDirs = allDirs([allDirs.isdir]);

    % Filter to stack_* dirs (skip '.' and '..')
    stackDirs = {};
    for i = 1:numel(allDirs)
        name = allDirs(i).name;
        if strcmp(name, '.') || strcmp(name, '..')
            continue;
        end
        if startsWith(name, 'stack_')
            stackDirs{end+1} = fullfile(allDirs(i).folder, name); %#ok<AGROW>
        end
    end

    fprintf('Found %d stack folders under:\n  %s\n', numel(stackDirs), masterDir);

    % Regex: prefix (Cam_*_), timepoint (5 digits), suffix (extension, including .lux.h5 etc)
    expr = '^(Cam_[^_]+_)(\d{5})(\..+)$';

    for d = 1:numel(stackDirs)
        stackDir = stackDirs{d};
        fprintf('\nProcessing: %s\n', stackDir);

        listing = dir(stackDir);
        listing = listing(~[listing.isdir]); % files only

        % Collect matching files and their old timepoints
        matchFiles = {};
        oldTps = [];

        for k = 1:numel(listing)
            fname = listing(k).name;
            tok = regexp(fname, expr, 'tokens', 'once');
            if ~isempty(tok)
                matchFiles{end+1} = fname; %#ok<AGROW>
                oldTps(end+1) = str2double(tok{2}); %#ok<AGROW>
            end
        end

        if isempty(matchFiles)
            fprintf('  (No matching Cam_*_#####.* files)\n');
            continue;
        end

        % Sort by old timepoint descending to avoid collisions
        [~, order] = sort(oldTps, 'descend');
        matchFiles = matchFiles(order);

        for k = 1:numel(matchFiles)
            oldName = matchFiles{k};
            tok = regexp(oldName, expr, 'tokens', 'once');
            prefix = tok{1};
            oldTp  = str2double(tok{2});
            suffix = tok{3};

            newTp = oldTp + desiredFirstTimepoint;
            newName = sprintf('%s%05d%s', prefix, newTp, suffix);

            oldPath = fullfile(stackDir, oldName);
            newPath = fullfile(stackDir, newName);

            if exist(newPath, 'file')
                error('Refusing to overwrite existing file:\n  %s\nwhile renaming:\n  %s', newPath, oldPath);
            end

            movefile(oldPath, newPath);
            fprintf('  %s -> %s\n', oldName, newName);
        end
    end

    fprintf('\nDone.\n');
end