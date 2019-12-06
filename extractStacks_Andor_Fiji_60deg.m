%% script to disentangle micro manager data; 
addpath(genpath('/home/streicha/code/imsane/imsane'))
 
 addpath('~/code/spimCode/')
nFiles     = 1; % automatize determination of number of tif files in folder;
nPos = 3;  
if nFiles <=1
     d = dir;
     temp = (length(d)-2)/nPos;
     if (temp-round(temp)) == 0
         nFiles = temp;    
     end
end
 
 
pos = 0:(nPos-1);
 
%nFilesMeta = min(nFiles,5); % number of files to read to determine metadata from;
%nameDummyPos  = 'TimeBeads_1_MMStack_Pos';
flag = 0;
for k = 1 : length(d)
    if (length(d(k).name) > 9)
        if strcmp(d(k).name(end-7:end),'.ome.tif') 
            flag = 1; 
            nameDummyPos = d(k).name(1:end-9);
            break;
        end
    end
end
if flag == 0
nameDummyPos  = 'TimeRot_1_MMStack_Pos';
end
% determine the meta data of the ome files; 
p = 1;
nameDummy = [nameDummyPos,num2str(pos(p))];
%meta = getMetaData(nameDummy,nFilesMeta);
meta(1).stacksize = 173;
meta(1).nchan = 1;
save('Meta.mat','meta');
%%
name  = [nameDummy,'.ome.tif'];
image = imread(name,1); 
% currently we assume two cameras are used;
%image = image(1:2:end,1:2:end);
image = zeros(size(image,1),size(image,2),meta(1).stacksize,class(image));
image2 = image;
%
counter   = 1;
col       = 1;
ncol      = meta(1).nchan;
TimeCount = 0;
 
mkdir('data');
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
%
for p = 1 : length(pos)
    disp(['Running pos ',num2str(pos(p))]);
    % loop over positiosn; the first one will be stack00, aso.
    nameDummy = [nameDummyPos,num2str(pos(p))];
    counter   = 1;
col       = 1;
TimeCount = 0;
    for time = TInit:TEnd
        disp(['Current time is ',num2str(time), ' of total ', num2str(TEnd)]);    
        if time > 0 
            name = [nameDummy,'_',num2str(time),'.ome.tif'];
        else 
            name = [nameDummy,'.ome.tif'];
        end
        
        filesize = length(imfinfo(name));
        for k = 1 : filesize%/2;      
            temp =  imread(name,k )- 100;
            image(:,:,counter)  = temp(1:end,1:end);
 
            counter = counter+1;
            if counter > 2*meta(1).stacksize
               image2 = image(:,end:-1:1,2:2:end);
               image2 = image2(:,:,end:-1:1);% for fusion using fiji image2 stack must be reversed in z. 
               
               %image2 = image2(end:-1:1,:,:); % when the image is flipped
               image2 = image2(:,end:-1:1,:);
               image  = image(:,:,1:2:end);
               counter = 1;
 
               %name1  = ['data\LC\stack',num2str(p-1),'\Time_',sprintf('%06d',TimeCount),'_c',num2str(col),'.tif'];
               %name2  = ['data\RC\stack',num2str(p-1),'\Time_',sprintf('%06d',TimeCount),'_c',num2str(col),'.tif'];
               
               name1 = ['data/Time_',sprintf('%06d',TimeCount),'_Angle_',num2str((p-1)*60),'_c',num2str(col),'.ome.tif'];
               name2 = ['data/Time_',sprintf('%06d',TimeCount),'_Angle_',num2str((p-1)*60+180),'_c',num2str(col),'.ome.tif'];
               %for j = 1 : meta(1).stacksize
               %    imwrite( image(:,:,j)',name1,'Compression','none','WriteMode','Append');
               %    imwrite(image2(:,:,j)',name2,'Compression','none','WriteMode','Append');
               %end
               % much faster saving than imwrite
               
               %disp('pause');
               
               %saveastiff(image,name1,options);
               %saveastiff(image2,name2,options);
              %  im = image(1:4:end,1:4:end,:);
               im = image;
               im = reshape(im,[size(im,1),size(im,2),1,size(im,3),1]); 
               bfsave(im,name1,'dimensionOrder', 'XYTZC')
              % im = image2(1:4:end,1:4:end,:);
               im = image2;
               im = reshape(im,[size(im,1),size(im,2),1,size(im,3),1]);
               bfsave(im,name2,'dimensionOrder', 'XYTZC')
               
               col = col+1;
               if col > ncol %&& p == length(pos)
                   TimeCount = TimeCount+1;
                   col = 1;
               end
 
 
            end
        end
        
        disp(['    Done with time ',num2str(time)]);
    end
    
     disp(['Done with pos ',num2str(pos(p))]);
end
 
disp('done')
