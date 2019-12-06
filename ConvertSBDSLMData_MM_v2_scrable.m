addpath(genpath('~/code/spimCode'));
nFiles     = 1; % automatize determination of number of tif files in folder;


nFilesMeta = 5; % number of files to read to determine metadata from;
nameDummy  = 'goodembryo_crop.tif';'TimeCol_1_MMStack';

% determine the meta data of the ome files; 
%meta = getMetaData(nameDummy,nFilesMeta);
%save('Meta.mat','meta');
name  = 'goodembryo_crop.tif';%[nameDummy,'.ome.tif'];
image = imread(name,1);
% currently we assume two cameras are used;
meta(1).stacksize = 318;
meta(1).nchan = 1;
image = zeros(size(image,1),size(image,2),meta(1).stacksize,class(image));
image2 = image;

counter   = 1;
col       = 1;
ncol      = meta(1).nchan;
TimeCount = 0;

mkdir('fused');
mkdir('thumbs');
mkdir('thumbs/x');
mkdir('thumbs/y');
mkdir('thumbs/z');
TInit = 0;
TEnd = nFiles-1;
%%
for time = TInit:TEnd
    
    if time > 0 
        name = [nameDummy,'_',num2str(time),'.ome.tif'];
    end
    
    filesize = length(imfinfo(name));
    for k = 1 : filesize%/2;
        
        
        %image(:,:,counter)  = imread(name, 2*(k-1)+1 );
        image(:,:,counter)  = imread(name,k )- 100;
        % need to also read the second image; 
        
        %temp = imread(name,2*k);
        %image2(:,:,counter) = temp(:,end:-1:1); % flip the image;
        
        %disp('Not yet reading second camera image!');
        counter = counter+1;
        if counter > 2*meta(1).stacksize
            
           image2 = image(end:-1:1,:,2:2:end);
           image  = image(:,:,1:2:end);
           counter = 1;
           if col == 1 %time == TInit && col == 1
               disp('first time');
               %pause
               leftflag = 1;
               [weights,seg] = WeightImage(image,leftflag);
               
               leftflag = 0;
               [weights2,seg2] = WeightImage(image2,leftflag);
               
               % compute image weights before fusion
           end
           % first we should fuse the two different views to a single view
           % here;
           % write to file, and potentially increase the color counter;
           dynoff  = (median(image(seg>0)))/2+(2*min(image(seg>0)))/4;
           dynoff2 = (median(image2(seg2>0)))/2+(2*min(image2(seg2>0)))/4;
           
           % fuse the two images to obtain a single high quality image; 
           
           %if TimeCount == 2
           %    disp('stop it!');
           %    pause
           %end
           
           image = FuseImages(image-dynoff,image2-dynoff2,weights,weights2);
           
           name1  =  ['fused/Time_',sprintf('%06d',TimeCount),'.tif'];
           for j = 1 : meta(1).stacksize
               imwrite(image(:,:,j)',[name1(1:end-4),'_f_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
               %imwrite(image2(:,:,j),[name1(1:end-4),'_right_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
           end
           imwrite(max(image,[],3)',['thumbs/z/',name1(7:end-4),'_max_f_c',num2str(col),'_maxz.tif'],'Compression','none','WriteMode','Append');
           imwrite(squeeze(max(image,[],1)),['thumbs/x/',name1(7:end-4),'_max_f_c',num2str(col),'_maxx.tif'],'Compression','none','WriteMode','Append');
           imwrite(squeeze(max(image,[],2)),['thumbs/y/',name1(7:end-4),'_max_f_c',num2str(col),'_maxy.tif'],'Compression','none','WriteMode','Append');
           %imwrite(max(image2,[],3),[name1(1:end-4),'_max_right_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
           %imwrite(squeeze(max(image2,[],1)),['thumbs\',name1(7:end-4),'_max_right_crop_c',num2str(col),'_maxx.tif'],'Compression','none','WriteMode','Append');
           %imwrite(squeeze(max(image2,[],2)),['thumbs\',name1(7:end-4),'_max_right_crop_c',num2str(col),'_maxy.tif'],'Compression','none','WriteMode','Append');
                          
           col = col+1;
           if col > ncol
               TimeCount = TimeCount+1;
               col = 1;
           end
            
            
        end
    end
    time
end