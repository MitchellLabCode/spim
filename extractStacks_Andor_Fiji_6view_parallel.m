%% script to disentangle micro manager data; 


%!!!!important!!!! 
% - always check that nFiles corresponds to the number of data
%   files per position in the raw data, unless set to 1 for automated
%   determination.
% - check the parameters, including angle_increment, bin, and ls

% To automatize determination of number of tif files in folder to read 
% (which is different than number to write, in general), set nFiles = 1.
nFiles     = 1; 
pos = [0,1,2,3,4,5];  % List of camera views
angle_increment = 30 ;  % How far to advance in degrees from view to view

if nFiles <=1
    % npm: We want to count # .ome files here, changing dir arg
    % and removing the '- 2' from '(length(d) - 2) / numel(pos)'
     d = dir('*.ome*');
     temp = (length(d))/numel(pos);
     if (temp-round(temp)) == 0
         nFiles = temp;  
     else 
         nFiles = floor(temp);
     end
end

nFilesMeta = min(nFiles,1); % number of files to read to determine metadata from;
flag = 0;
for k = 1 : length(d)
    if (length(d(k).name) > 9)
        if strcmp(d(k).name(end-7:end),'.ome.tif') 
            flag = 1; 
            % nameDummyPos might be something like 'TimeBeads_1_MMStack_Pos'
            nameDummyPos = d(k).name(1:end-9);
            break;
        end
    end
end
if flag == 0
nameDummyPos  = 'Sample2_3D_1_MMStack_Pos';
end
% determine the meta data of the ome files; 
p = 1;
nameDummy = [nameDummyPos,num2str(pos(p))];
meta = getMetaData(nameDummy,nFilesMeta);
save('Meta.mat','meta');

%%
name  = [nameDummy,'.ome.tif'];
imageDummy = imread(name,1); 
% Get stacksize for image preallocation
meta_stacksize = meta(1).stacksize ;

% Parameters (check before running)
ls = 1;  % the number of lightsheets
bin = 1;  % The binning amount in data acquisition
counter   = 1;
col       = 1;
ncol      = 1; % meta(1).nchan; %Set number of colors
TimeCount = 0;
datafolder = 'data';
%log file
fid = fopen('log.txt','a');
fprintf(fid,'%s\r\n',datetime);
fprintf(fid,'-------------------\r\n');
fprintf(fid,'%s\r\n',datafolder);
fprintf(fid,'number of colors: %i \r\n',ncol);
fprintf(fid,'number of positions: %i \r\n',numel(pos));
fprintf(fid,'binnig: %i\r\n',bin);
fprintf(fid,'\r\n\r\n');
fclose(fid);

Channel_Tag = 'Multi Camera-CameraChannelName';

mkdir(datafolder);
%mkdir('data\LC');
%mkdir('data\LC\stack0');
%mkdir('data\LC\stack1');
%mkdir('data\RC');
%mkdir('data\RC\stack0');
%mkdir('data\RC\stack1');
TInit = 0;
TEnd = nFiles-1;
%
clear options
options.append = true;
options.compress = 'no';

% Iterate over each camera angle
parfor p = 1: length(pos)
    % Create an image which is 3D and has same xy resolution as imread image
    image = zeros(size(imageDummy,1),size(imageDummy,2),meta_stacksize,class(imageDummy));
    image2 = image;

    disp(['Running pos ',num2str(pos(p))]);
    % loop over positiosn; the first one will be stack00, aso.
    nameDummy = [nameDummyPos,num2str(p-1)];
    counter   = 1;
    counter_1 = 1;
    counter_2 = 1; %% we need to count for both cameras now
    col       = 1;
    TimeCount = 0;
    for time = TInit:TEnd
        disp(['Current time is ',num2str(time), ' of total ', num2str(TEnd)]);    
        if time > 0 
            name = [nameDummy,'_',num2str(time),'.ome.tif'];
        else 
            name = [nameDummy,'.ome.tif'];
        end
        
        % Preallocate arrays for which image is from which camera (cam)
        bla = imfinfo(name);
        cam = zeros(size(bla));
        filesize = length(bla);
        % determine the camera for each file; 
        for k = 1 : filesize
            val = bla(k).UnknownTags(end).Value;
            ind = strfind(val,Channel_Tag);
            key = val((ind+length(Channel_Tag)+3):(ind+length(Channel_Tag))+30);
            if ~contains(key,'-left')
               cam(k) = 1;
            else
               cam(k) = 2;
            end
        end
                
        for k = 1 : filesize     
            temp =  imread(name,k )- 100; % read image, and take away constant offset
            if cam(k) == 1 % which camera does the image come from
                image(:,:,counter_1)  = temp(1:end,1:end);
                counter_1 = counter_1+1;
            elseif cam(k) ==2 
                image2(:,:,counter_2) = temp(1:end,1:end);
                counter_2 = counter_2+1;
            end
            counter = counter+1;
            % Once we read 2 stacks, we write to disk for this timestamp
            if counter > 2* meta_stacksize
                image2 = image2(:,:,end:-1:1); % for fusion using fiji image2 stack must be reversed in z. 
                counter = 1;
                counter_1 = 1;
                counter_2 = 1;
                name1 = [datafolder,'\Time_',sprintf('%06d',TimeCount),'_Angle_',num2str((pos(p))*angle_increment+0*180),'_c',num2str(col),'_ls_',num2str(ls),'.ome.tif'];
                name2 = [datafolder,'\Time_',sprintf('%06d',TimeCount),'_Angle_',num2str((pos(p))*angle_increment+1*180),'_c',num2str(col),'_ls_',num2str(ls),'.ome.tif'];

                % Note that saving using imwrite (below) is slow, and 
                % saveastiff is also not as optimal, so instead use
                % Bio-Formats library and bsave
                
                % SAVE CAMERA 1
                % Handle case where data was binned
                if bin == 2
                    im = uint16( ( double( image(1:2:end,1:2:end,:))+double(image(2:2:end,1:2:end,:))+...
                         double(image(1:2:end,2:2:end,:))+double(image(2:2:end,2:2:end,:) ) )/4);
                elseif bin == 8
                  im=zeros(size(image,1)/8,size(image,2)/8,size(image,3));
                  for ii = 1:8
                      for jj=1:8
                        im = im + double(image(ii:8:end,jj:8:end,:));
                      end
                  end

                  im =uint16(im/64);
                else
                  im = image; 
                end %im = image(:,:,end:-1:1);
               im = reshape(im,[size(im,1),size(im,2),1,size(im,3),1]); 
               bfsave(im, name1, 'dimensionOrder', 'XYTZC')
               
               % SAVE CAMERA 2
               % Now do the second image, which is the second camera - npm
               if bin == 2
                 im =uint16( ( double( image2(1:2:end,1:2:end,:))+double(image2(2:2:end,1:2:end,:))+...
                     double(image2(1:2:end,2:2:end,:))+double(image2(2:2:end,2:2:end,:) ) )/4); 
               elseif bin == 8
                 im=zeros(size(image2,1)/8,size(image2,2)/8,size(image2,3));
                 for ii = 1:8
                    for jj=1:8
                      im = im + double(image2(ii:8:end,jj:8:end,:));
                    end
                 end 
                 im = uint16(im/64);
               else 
                 im = image2; 
               end 
               
               im = reshape(im,[size(im,1),size(im,2),1,size(im,3),1]);
               % bfsave saves a 5D matrix into OME-TIFF using Bio-Formats
               bfsave(im, name2, 'dimensionOrder', 'XYTZC')
               
               col = col+1;
               % Check if we have filled up this stack (all columns)
               if col > ncol 
                    % This stack is done, so start a new timestamp and
                    % return to column 1
                    disp([' -> Done with output image timestamp ', num2str(TimeCount)]);

                    % Add to log: append to log.txt using 'a' option
                    fid = fopen('log.txt','a');
                    fprintf(fid, ['  -> Done with output image timestamp ', num2str(TimeCount)]);
                    fprintf(fid, '\r\n');
                    fclose(fid);

                    TimeCount = TimeCount+1;
                    col = 1;
               end
            end
        end
        disp(['    Done with file timestamp ',num2str(time)]);
        
        % Add to log: append to log.txt using 'a' option
        fid = fopen('log.txt','a');
        fprintf(fid, ['    Done with filetimestamp ', num2str(time)]);
        fprintf(fid, '\r\n');
        fclose(fid);
    end
    disp(['Done with pos ', num2str(pos(p))]);
    
    % Add to log: append to log.txt using 'a' option
    fid = fopen('log.txt','a');
    fprintf(fid, ['Done with pos ', num2str(pos(p))]);
    fprintf(fid, '\r\n');
    fclose(fid);
end

disp('done')
