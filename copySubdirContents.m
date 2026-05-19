function copySubdirContents(sourceDir, destDir)
    % copySubdirContents Copies contents of subdirectories from source to destination.
    %
    % Inputs:
    %   sourceDir - String or char array of the main source directory path
    %   destDir   - String or char array of the main destination directory path

    % Verify that both main directories exist
    if ~isfolder(sourceDir)
        error('Source directory does not exist: %s', sourceDir);
    end
    if ~isfolder(destDir)
        error('Destination directory does not exist: %s', destDir);
    end

    % Get a list of all items in the source directory
    items = dir(sourceDir);

    % Filter out non-directories and the standard '.' and '..' folders
    isDir = [items.isdir];
    dirNames = {items.name};
    validDirs = items(isDir & ~ismember(dirNames, {'.', '..'}));

    % Loop through each valid subdirectory
    for i = 1:length(validDirs)
        subDirName = validDirs(i).name;
        
        % Construct full paths for the source and destination subdirectories
        sourceSubDir = fullfile(sourceDir, subDirName);
        destSubDir   = fullfile(destDir, subDirName);
        
        % Check if the matching subdirectory exists in the destination
        if isfolder(destSubDir)
            fprintf('Copying contents from: %s\n                 To: %s\n', subDirName, destSubDir);
            
            % Copy all contents (*) of the source subdirectory to the destination
            % The 'f' flag forces the copy (overwriting read-only files if necessary)
            [success, msg, msgID] = copyfile(fullfile(sourceSubDir, '*'), destSubDir, 'f');
            
            if ~success
                warning('Failed to copy contents for %s. Reason: %s', subDirName, msg);
            end
        else
            % If you prefer to create the missing directories automatically, 
            % you can replace the warning below with:
            % mkdir(destSubDir); copyfile(fullfile(sourceSubDir, '*'), destSubDir, 'f');
            
            fprintf('Skipping %s: Matching subdirectory not found in destination.\n', subDirName);
        end
    end
    
    fprintf('Copy operation complete.\n');
end