
nFiles     = 70; % automatize determination of number of tif files in folder;


nFilesMeta = 5; % number of files to read to determine metadata from;
nameDummy  = 'TimeCol_1_MMStack';

% determine the meta data of the ome files; 
meta = getMetaData(nameDummy,nFilesMeta);

name  = [nameDummy,'.ome.tif'];
image = imread(name,1);
% currently we assume two cameras are used;
image = zeros(size(image,1),size(image,2),meta(1).stacksize,class(image));
image2 = image;

counter   = 1;
col       = 1;
ncol      = meta(1).nchan;
TimeCount = 0;

mkdir('fused');
mkdir('thumbs');

%%
for time = 0 : nFiles-1
    
    if time > 0 
        name = [nameDummy,'_',num2str(time),'.ome.tif']
    end
    
    filesize = length(imfinfo(name));
    
    for k = 1 : filesize/2;
        
        image(:,:,counter)  = imread(name,2*(k-1)+1);
        % need to also read the second image; 
        %image2(:,:,counter) = imread(name,2*k);
        disp('Not yet reading second camera image!');
        counter = counter+1;
        if counter > meta(1).stacksize
            
            counter = 1;
            % write to file, and potentially increase the color counter;
            
           name1  =  ['fused\Time_',sprintf('%06d',TimeCount),'.tif'];
           for j = 1 : meta(1).stacksize
               imwrite(image(:,:,j),[name1(1:end-4),'_left_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
               %imwrite(image2(:,:,j),[name1(1:end-4),'_right_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
           end
           imwrite(max(image,[],3),['thumbs\',name1(7:end-4),'_max_left_crop_c',num2str(col),'_maxz.tif'],'Compression','none','WriteMode','Append');
           imwrite(squeeze(max(image,[],1)),['thumbs\',name1(7:end-4),'_max_left_crop_c',num2str(col),'_maxx.tif'],'Compression','none','WriteMode','Append');
           imwrite(squeeze(max(image,[],2)),['thumbs\',name1(7:end-4),'_max_left_crop_c',num2str(col),'_maxy.tif'],'Compression','none','WriteMode','Append');
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