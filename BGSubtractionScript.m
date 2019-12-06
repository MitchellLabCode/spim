addpath(genpath('/home/streicha/code'))
nameDummy = 'Time0.tif';%'TimeCol_7_MMStack.ome.tif';
%cd fused
mkdir mips
BallRad   = 40;
IBd       = [0, 2^16-1];
aspect    = .26/.26;
%
for t =  0 
    
    name = sprintf(nameDummy,t);
    
    stackSize = length(imfinfo(name));
    image = imread(name,1);
    image = image(1:1:end,1:1:end);
    kInd = 1 : 1 : stackSize;
    for k = 1 : length(kInd)
        
        temp = imread(name,kInd(k));
        temp = medfilt2(temp,[4 4]);
        temp = temp(1:1:end,1:1:end);
        image(:,:,k) = temp-imopen(temp,strel('disk',14))/2;
    end
    %image = image(1:2:end,1:2:end,:);
    tic
    imagebgc = TopHatFFT(image,aspect,BallRad,IBd);
    toc
    
    level = graythresh(mat2gray(imagebgc));
    seg = mat2gray(imagebgc)>level/2;
    
    seg = bwareaopen(seg,20);
    
    imagebgc = imagebgc.*uint16(seg);
    
    %for k = 1 : stackSize
        
    %    imwrite(imagebgc(:,:,k),[name(1:end-4) '_bgc.tif'],'Compression','none','WriteMode','append');
    %end
    %imwrite(max(imagebgc,[],3),[name(1:end-4) '_bgc_mip.tif'],'Compression','none','WriteMode','append');
    imwrite(max(imagebgc(:,:,1:round(end/2)),[],3),['mips/',name(1:end-4) '_bgc_mip1.tif'],'Compression','none','WriteMode','overwrite');
    imwrite(max(imagebgc(:,:,round(end/2):end),[],3),['mips/',name(1:end-4) '_bgc_mip2.tif'],'Compression','none','WriteMode','overwrite');
    imwrite(squeeze(max(imagebgc(1:round(end/2),:,:),[],1)),['mips/',name(1:end-4) '_bgc_mip_y1.tif'],'Compression','none','WriteMode','overwrite');
    imwrite(squeeze(max(imagebgc(round(end/2):end,:,:),[],1)),['mips/',name(1:end-4) '_bgc_mip_y2.tif'],'Compression','none','WriteMode','overwrite');
end