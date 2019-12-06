datadir = '/mnt/crunch/how24bGal4UasHistRFP/201903181645_good_lotsofbeads/Time6views_60sec_1p4um_25x_2p0mW_exp0p50_6' ;
searchstr = 'data*Time_0*_Angle_*.ome.tif' ;
removestr = 'data\' ;
replacementstr = '' ;
% Rename all files matching string pattern inside a directory. 
% Note that this can only replace the filename, not other parts of the path
% or the file extension.
% npm 2019

% Retrieve all the matching file names
files = dir(fullfile(datadir, searchstr)) 

for i=1:length(files)
    % split up the filename
    % [pathn, filename, extension] = fileparts(files(i).name) 
    pathn = files(i).folder ;
    filename = files(i).name ;
    
    % Create the new name
    newfilename = strrep(filename, removestr, replacementstr) ;
    oldfn = fullfile(pathn, filename) ;
    newfn = fullfile(pathn, newfilename) ;
    
    sprintf(['Replacing ', filename, '\n  with -> ', newfilename])
    % sprintf(['Replacing ', strrep(oldfn, '\', '\\'), 
    %          '  with -> ', strrep(newfn, '\', '\\')])
    % fprintf(['Replacing ', oldfn, '\n  with -> ', newfn, '\n'])
    
    movefile(oldfn, newfn)
end


return