function batch_rename_files(datadir, searchstr, removestr, replacementstr)
% Rename all files matching string pattern inside a directory. 
% Note that this can only replace the filename, not other parts of the path
% or the file extension.
% npm 2019

% Retrieve all the matching file names
files = dir(fullfile(datadir, searchstr))

for i=1:length(files)
    % split up the filename
    [pathn, filename, extension] = fileparts(files(i).name)
    % Create the new name
    newfn = strrep(filename, removestr, replacementstr)
    movefile([pathn filename extension], [pathn newfn extension])
end


return