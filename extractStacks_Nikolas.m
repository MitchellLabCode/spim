%% script to disentangle micro manager data; 
% Micro manager saves data as fast as possible into a series of .ome.tif
% files, which contain as many pages as allowed by the filesystem (4GB).
% Each page corresponds to a certain time, z-coordinate, and camera.
% Pages are ordered acccording to time and z-coordinate, but we need to
% figure out where each timepoint (i.e. each full z-stack) begins and ends,
% and disentangle the pages taken by cameras 1 and 2 at the same time.

% log4j warnings can be ignored - they are due to the bioformats library
% used for saving images. bioformats' bfsave is much faster than the matlab
% alternatives imwrite and saveastiff.

% Output files have the following name structure:
% Time_{tttttt}_Angle_{a}_c_{c}_ls_{ls}_.ome.tif
% where 't' is the time, 'a' is the view angle and 'c' is the color channel

addpath(genpath('/mnt/data/code/spimCode'));
addpath('/mnt/data/code/imsaneV1/external/bfmatlab');

%% define parameters

%!!!!important!!!! : check the parameters

pos = [0, 1, 2, 3]; % list of camera views.
angle_increment = 45; % angle increment for output filename
datafolder = 'data2'; % folder for saving results
ncol = 1;  % Set number of colors
bin = 1; % Bin images before saving (i.e. downsample by summing together
% adjacent pixels) if bin = 1, no binning is performed. In case bin = 2, 8,
% etc., saved images are downsampled by a factor of 4 resp. 64.
ls = 1; % random parameter for output filenames
background_intensity = 100; % subtract background intensity. For some
% reason, images come out with a minimal intensity of 100 instead of 0.

%% determine number of files and read metadata

% determine number of files.
file_list = dir('*.ome.tif');
nFiles = floor(length(file_list) / length(pos));

% then determine the meta data of the ome files by reading data from one 
% arbitrary file (should be the same for all files)
filename_prefix = file_list(1).name(1:end-9);
meta = getMetaData([filename_prefix, num2str(pos(1))], 1);
save('Meta.mat','meta');

%% Main loop

% We assume two cameras are used -> each position contains two images.

% Create 2 empty arrays to hold the results
image1 = imread([filename_prefix, num2str(pos(1)), '.ome.tif'], 1); 
image1 = zeros(size(image1,1), size(image1,2), meta(1).stacksize, class(image1));
image2 = image1;

% write parameters to log file
fid = fopen('log.txt','a'); % log file
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

% iterate over positions;
for p = pos
    disp(['Running pos ', num2str(p)]);
    
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
    
    % z_counter etc determine the z-coordinate of a page. start at 1,
    % since they matlab indices.
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
    for file_number = 0 : nFiles-1 % used to be -2
        disp(['Currently processing file ',num2str(file_number), ' of total ', num2str(nFiles-1)]);
        % get name of file
        if file_number > 0 
            name = [filename_prefix, num2str(p), '_',num2str(file_number),'.ome.tif'];
        else 
            name = [filename_prefix, num2str(p), '.ome.tif'];
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
        for k = 1 : number_pages
            % read image (removed: "-100")
            if cam(k) == 1 % which camera does the image come from
                image1(:, :, z_counter_1) = imread(name, k)-background_intensity;
                z_counter_1 = z_counter_1+1;
            elseif cam(k) == 2 
                image2(:,:,z_counter_2) = imread(name, k)-background_intensity;
                z_counter_2 = z_counter_2+1;
            end
            z_counter = z_counter+1;
            
            % Instead of assuming equipartition of camera 1 and camera 2
            % within a 2*stackSize stack, we have to check each camera
            % individually. (Bug Identified in Jan 2020)
            
            % If we have read a complete z-stack (one complete timepoint),
            % save it and reset the z-counters. Do this for cameras
            % 1 and 2 separately
            if z_counter_1 > meta(1).stacksize % that's the number of z-stacks
                name1 = [datafolder,'/Time_',sprintf('%06d',t_counter_1),'_Angle_',...
                    num2str(p*angle_increment+0*180),'_c',num2str(col),'_ls_',num2str(ls),'.ome.tif'];
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
                bfsave(image_to_save, name1, 'dimensionOrder', 'XYTZC')
            end
            if z_counter_2 > meta(1).stacksize
                name2 = [datafolder,'/Time_',sprintf('%06d',t_counter_2),'_Angle_',...
                    num2str(p*angle_increment+1*180),'_c',num2str(col),'_ls_',num2str(ls),'.ome.tif'];
                z_counter_2 = 1;
                t_counter_2 = t_counter_2+1;
                if bin > 1
                    image_to_save = binning(image2, bin);
                else
                    image_to_save = image2;
                end
                image_to_save = image_to_save(:,:,end:-1:1);
                % for fusion using fiji, image2 stack must be reversed in z. 
                image_to_save = reshape(image_to_save,[size(image_to_save,1),...
                    size(image_to_save,2),1,size(image_to_save,3),1]);
                % above line was missing, creating incorrect metadata for
                % all images taken from camera 2. Identified 09/14/21.
                bfsave(image_to_save, name2, 'dimensionOrder', 'XYTZC')
            end
            
            % now think about colour; MAKE SURE IS CORRECT!
            if z_counter > 2*meta(1).stacksize
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
    end
    disp(['Done with pos ', num2str(p)]);
    fid = fopen('log.txt','a');
    fprintf(fid, ['Done with pos ', num2str(p)]);
    fclose(fid);
end


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

    name1 = [datafolder,'/Time_',sprintf('%06d',t_counter),'_Angle_',num2str((p-1)*angle_increment+0*180),'_c',...
        num2str(col),'_ls_',num2str(ls),'.ome.tif'];
    name2 = [datafolder,'/Time_',sprintf('%06d',t_counter),'_Angle_',num2str((p-1)*angle_increment+1*180),'_c',...
        num2str(col),'_ls_',num2str(ls),'.ome.tif'];

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
