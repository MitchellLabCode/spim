%% script to disentangle micro manager data and crop out bad z-slices; 
% Micro manager saves data as fast as possible into a series of .ome.tif
% files, which contain as many pages as allowed by the filesystem (4GB).
% Each page corresponds to a certain time, z-coordinate, and camera.
% Pages are ordered acccording to time and z-coordinate, but we need to
% figure out where each timepoint (i.e. each full z-stack) begins and ends,
% and disentangle the pages taken by cameras 1 and 2 at the same time.

% log4j warnigit_repo_clone/imsane/imsane/externalngs can be ignored - they are due to the bioformats library
% used for saving images. bioformats' bfsave is much faster than the matlab
% alternatives imwrite and saveastiff.

% Output files have the following name structure:
% Time_{tttttt}_Angle_{a}_c_{c}_ls_{ls}_.ome.tif
% where 't' is the time, 'a' is the view angle and 'c' is the color channel

% This version can be run in parallel to imaging, because it automatically
% waits until a .ome.tif file is complete.
% If desired, we can also crop of n z slices at the beginning and end
% of each image since they can be corrputed.

close all
clearvars

% addpath('C:\Users\lichtblatt\Desktop\bfmatlab');
% addpath('C:\Users\lichtblatt\Desktop\SPIMcode-main\SPIMcode-main');
addpath('/mnt/data/code/SPIMcode')
addpath('/mnt/data/code/tubular/utility/bfmatlab/')

%% define parameters

% crop away unnecessary z-slices
% to do: cjppse n_zcrop as a function of camera
n_zcrop = [0 0];

%!!!!important!!!! : check the parameters

manual_n_z = 482; % number of z stacks. if -1, tries automatic detection
% from the metadata, this fails for larger (>255) z pages

pos = [0,1,2,3]; % 1, 2, 3]; % list of dual-camera angles.
angle_increment = 45; % angle increment for output filename
datafolder = '/mnt/crunch/actb2-mem-cherry_h2afva-gfp/202307111346_4views_1p4um_1ms/'; % folder for saving results
ncol = 2;  % Set number of colors
bin = 1; % Bin images before saving (i.e. downsample by summing together
% adjacent pixels) if bin = 1, no binning is performed. In case bin = 2, 8,
% etc., saved images are downsampled by a factor of 4 resp. 64.
ls = 1; % random parameter for output filenames
background_intensity = 0; % subtract background intensity. For some
% reason, images come out with a minimal intensity of 100 instead of 0.

%% determine number of files and read metadata

% determine number of files. You can define this manually
% to the expected number of files if you want to run it in parallel
% leave at 0 for automatic detection, or at -1 if you want 
% to run until wait timeout happens

file_list = dir('*.ome.tif');
nFiles = 0;
if nFiles == 0
    nFiles = floor(length(file_list) / length(pos));
end
if nFiles == -1
    nFiles = 1000; % arbitrary, large number
end

max_wait = 3; % maximum time you are ok waiting for next .ome.tif, in minutes
t_pause = 3; % time you wait for each new .ome.tif, not important.
% script terminates after waiting max_wait

% if run live, the script waits n seconds for the next file
% to be written. if it does mot happen, finish.

% then determine the meta data of the ome files by reading data from one 
% arbitrary file (should be the same for all files)
filename_prefix = file_list(1).name(1:end-9);
meta = getMetaData([filename_prefix, num2str(pos(1))], 1);
if ~exist('Meta.mat', 'file')
    save('Meta.mat','meta');
else
    posstr = 'pos' ;
    for dmyi = 1:length(pos)
        posstr = [posstr num2str(pos(dmyi))] ;
    end
    save(['Meta_' posstr '.mat'])
end

%% Main loop

% We assume two cameras are used -> each position contains two images.

% write parameters to log file
if ~exist('log.txt', 'file')
    logfn = 'log.txt' ; % log file
else
    posstr = 'pos' ;
    for dmyi = 1:length(pos)
        posstr = [posstr num2str(pos(dmyi))] ;
    end
    logfn = ['log_' posstr '.txt'] ;
end
fid = fopen(logfn, 'a') ;
fprintf(fid,'-------------------\n');
fprintf(fid,'%s\n', datafolder);
fprintf(fid,'number of colors: %i\n', ncol);
fprintf(fid,'number of positions: %i\n', numel(pos));
fprintf(fid,'binnig: %i\n', bin);
fprintf(fid,'\n\n', bin);
fclose(fid);

Channel_Tag = 'Multi Camera-CameraChannelName';

mkdir(datafolder);
clear options
options.append = true;
options.compress = 'no';

%% parallel loop over positions

% iterate over positions;
%parpool(4)
for p = pos
    disp(['Running pos ', num2str(p) ...
        ': Angle=' num2str(p*angle_increment) ...
        ' & ' num2str(p*angle_increment+180)]);
    
    % sequentially read the pages written to each .ome.tif. Using the
    % counter variables, figure out what camera and timepoint a given page
    % belongs to. Pages are assumed to be ordered according to z-coord and
    % time, but not according to camera position. We have to figure out
    % where a complete z-stack (i.e. all the images corresponding to one
    % timepoint) begins and ends, and the save that. At the same time, we
    % need to separate the images taken by camera 1 and camera 2. 
    % The number of z-slices / time point is given by
    % meta(1).stacksize. Pages are sequentially read, and counter variables
    % used to determine whether a whole z-stack has been read. In this
    % case, the next page corresponds to the subsequent timepoint.
    % In this way, the pages are ordered and stored in the "image1" array.
    % Because there are two cameras, and the pages are not ordered by
    % camera in the .ome.tif files, there are two "image" array, and
    % two of each counter variable.

    if manual_n_z > -1
        n_z_stacks = manual_n_z;
    else
        n_z_stacks = meta(1).stacksize;
    end


    % Create 2 empty arrays to hold the results
    image1 = imread([filename_prefix, num2str(pos(1)), '.ome.tif'], 1); 
    image1 = zeros(size(image1,1), size(image1,2), n_z_stacks, class(image1));
    image2 = image1;

    % z_counter etc determine the z-coordinate of a page
    z_counter   = 1;
    z_counter_1 = 1;
    z_counter_2 = 1;
    % time_counter etc. determine the time point of a page.
    t_counter   = 0;
    t_counter_1 = 0;
    t_counter_2 = 0;
    % color counter - only important in the case of multiple channels /
    % timepoint
    col = 1;

    % loop over files
    wait_timeout = 0; % will be used to break out of loop in case of timeout
    for file_number = 0 : nFiles-1 % used to be -2
        disp(['Currently processing file ',num2str(file_number), ' of total ', num2str(nFiles-1)]);
        % get name of file
        if file_number > 0 
            name = [filename_prefix, num2str(p), '_',num2str(file_number),'.ome.tif'];
            name_nxt = [filename_prefix, num2str(p), '_',num2str(file_number+1),'.ome.tif'];
        else 
            name = [filename_prefix, num2str(p), '.ome.tif'];
            name_nxt = [filename_prefix, num2str(p), '_',num2str(1),'.ome.tif'];
        end
        % check if next file exists, so we know current .ome.tif complete. Else, wait a little
        file_exists = isfile(name_nxt) | (file_number == nFiles);
        % keep track of total wait time.
        t_wait_total = 0;
        while ~file_exists
            disp("waiting for file ...")
            pause(60*t_pause)
            t_wait_total = t_wait_total+t_pause;
            file_exists = isfile(name_nxt);
            if t_wait_total >= max_wait
                file_exists = 1;
                wait_timeout = 1;
            end
        end
        
        % for each page in file, determine the camera which recorded it
        info_struct = imfinfo(name);
        
        number_pages = length(info_struct);
        cam = zeros(number_pages);
        for k = 1 : number_pages
            val = info_struct(k).UnknownTags(end).Value;
            ind = strfind(val,Channel_Tag);
            key = val((ind+length(Channel_Tag)+3):(ind+length(Channel_Tag))+30);
            if ~contains(key, '-left')
               cam(k) = 1;
            else
               cam(k) = 2;
            end
        end
        
        % iterate over the pages in the .ome.tif file
        disp("reading ...")
        for k = 1 : number_pages
            % read image
            if cam(k) == 1 % which camera does the image come from
                image1(:, :, z_counter_1) = imread(name, k)-background_intensity;
                z_counter_1 = z_counter_1+1;
            elseif cam(k) == 2 
                image2(:,:,z_counter_2) = imread(name, k)-background_intensity;
                z_counter_2 = z_counter_2+1;
            end
            z_counter = z_counter+1;
            
            % If we have read a complete z-stack (one complete timepoint),
            % save it and reset the z-counters. Do this for cameras
            % 1 and 2 separately
            if z_counter_1 > n_z_stacks % that's the number of z-stacks
                name1 = fullfile(datafolder, ['Time_',sprintf('%06d',t_counter_1),'_Angle_',...
                    num2str(p*angle_increment+0*180),'_c',num2str(col),'_ls_',num2str(ls),'.ome.tif']);
                z_counter_1 = 1;
                t_counter_1 = t_counter_1+1;
                % optionally, use binning (=downsampling)
                if bin > 1
                    image_to_save = binning(image1, bin);
                else
                    image_to_save = image1;
                end
                % re-order dimensions so that we get a z-stack (and not a
                % t-stack) in fiji. Important for further steps in the 
                % pipeline to work properly.
                image_to_save = reshape(image_to_save,[size(image_to_save,1),...
                    size(image_to_save,2),1,size(image_to_save,3),1]);
                
                % if removing bad z-slices is desired, do it here:
                if sum(n_zcrop) > 0
                    image_to_save = image_to_save(:, :, :, 1+n_zcrop(1):end-n_zcrop(2));
                end
                disp("saving")
                bfsave(image_to_save, name1, 'dimensionOrder', 'XYTZC')
            end
            if z_counter_2 > n_z_stacks
                name2 = fullfile(datafolder, ['Time_',sprintf('%06d',t_counter_2),'_Angle_',...
                    num2str(p*angle_increment+1*180),'_c',num2str(col),'_ls_',num2str(ls),'.ome.tif']);
                z_counter_2 = 1;
                t_counter_2 = t_counter_2+1;
                if bin > 1
                    image_to_save = binning(image2);
                else
                    image_to_save = image2;
                end
                image_to_save = image_to_save(:,:,end:-1:1);
                % for fusion using fiji, image2 stack must be reversed in z. 
                image_to_save = reshape(image_to_save,[size(image_to_save,1),...
                    size(image_to_save,2),1,size(image_to_save,3),1]);
                % if removing bad z-slices is desired, do it here:
                if sum(n_zcrop) > 0
                    image_to_save = image_to_save(:,:,:,1+n_zcrop(1):end-n_zcrop(2));
                end
                disp("saving")
                bfsave(image_to_save, name2, 'dimensionOrder', 'XYTZC')
            end
            
            % now think about colour; MAKE SURE IS CORRECT!
            if z_counter > 2*n_z_stacks
                z_counter = 1;
                col = col+1;
                if col > ncol
                    t_counter = t_counter+1;
                    col = 1;
                end
            end
        end
        disp(['Done with file_number ',num2str(file_number)]);
        fid = fopen('log.txt','a');
        fprintf(fid, ['Done with file_number ',num2str(file_number)]);
        fclose(fid);
        if wait_timeout
            disp("Time out, assuming experiment finished")
            break
        end
    end
    disp(['Done with pos ', num2str(p)]);
    fid = fopen('log.txt','a');
    fprintf(fid, ['Done with pos ', num2str(p)]);
    fclose(fid);

    %% Fixing micro manager bug for last timepoint

    if z_counter ~=1 % meaning the last z-stack is not comletely filled
        disp('Micro manager missing file bug detected. Attempting fix.');
        fid = fopen('log.txt','a');
        fprintf(fid, 'Micro manager missing file bug detected. Attempting fix.');
        fclose(fid);
    
        % There is a potential bug with missing images from micro manager. 
        % This leads to partial missing files in the last timepoint of an
        % experiment. To overcome this, we force write the last potentially
        % incomplete stack.
    
        name1 = fullfile(datafolder,['Time_',sprintf('%06d',t_counter),'_Angle_',num2str((p-1)*angle_increment+0*180),'_c',...
            num2str(col),'_ls_',num2str(ls),'.ome.tif']);
        name2 = fullfile(datafolder, ['Time_',sprintf('%06d',t_counter),'_Angle_',num2str((p-1)*angle_increment+1*180),'_c',...
            num2str(col),'_ls_',num2str(ls),'.ome.tif']);
    
        % start with image1 (camera 1)
        if bin > 1
            image1 = binning(image1,bin);
        end
        image1 = reshape(image1,[size(image1,1),size(image1,2),1,size(image1,3),1]); 
        bfsave(image1,name1,'dimensionOrder', 'XYTZC')
        % now do the same thing for image2 (camera 2)
        image2 = image2(:,:,end:-1:1); % for fusion using fiji, image2 stack must be reversed in z. 
        if bin > 1
            image2 = binning(image2,bin);
        end
        image2 = reshape(image2,[size(image2,1),size(image2,2),1,size(image2,3),1]); 
        bfsave(image2,name2,'dimensionOrder', 'XYTZC')
    end
end

disp('Done')
fid = fopen('log.txt','a');
fprintf(fid, 'Done');
fclose(fid);


%% local function for binning

function img_binned = binning(img, bin)
% BINNING Local function that downsamples 3d-image by binning.
    img_binned = zeros(size(img,1)/bin,size(img,2)/bin,size(img,3));
    for ii = 1:bin
        for jj=1:bin
            img_binned = img_binned + double(img(ii:bin:end,jj:bin:end,:));
        end
    end
    img_binned = uint16(image_temp/bin^2);
end
