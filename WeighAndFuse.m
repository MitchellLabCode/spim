% Bead detection and point matching. 
% File and Folder Information
ExperimentFolderName = cd;%'G:\embryo\eve\20140606161200\TimeCol_1';

LFolderName = '';%'L/';
RFolderName = '';%'R/';

k = 1;
t = 0;

nch = 1;

Rotations = [0,pi/2]; % angles of different rotation views; 

% spatial resoution
res = [.26,.26,1.25]; %dx, dy, dz;
UntrackedPenalty = 100;
cd(ExperimentFolderName);
addpath(genpath('C:\Users\spim\Desktop\code\'))

%% 

if ~exist('fused','dir')
    mkdir('fused');
end
%% 
camOffset = 100; % General camera offset. 

for t = 37:119
    tic
  
        LFileName   = ['Time_',sprintf('%06d',t),'_left_crop_c1.tif'];
        RFileName   = ['Time_',sprintf('%06d',t),'_right_crop_c1.tif'];
        
        if nch == 2
        	LFileName2   = ['Time_',sprintf('%06d',t),'_left_crop_c1.tif'];
            RFileName2   = ['Time_',sprintf('%06d',t),'_right_crop_c1.tif'];
        else
            LFileName2 = LFileName;
            RFileName2 = RFileName;
        end
        
        StackSize = length(imfinfo([LFileName(1:end-3),'tif']));
        imageL    = imread([LFileName(1:end-3),'tif']);
        imageL    = zeros([size(imageL,1),size(imageL,2),StackSize],'uint16');
        imageR    = imageL;
        imageL2   = imageL;
        imageR2   = imageR;
        for z = 1 : size(imageL,3)
            imageL(:,:,z) = imread(LFileName,z)-camOffset;
            tmp = imread(RFileName,z)-camOffset;
            imageR(:,:,z) = tmp(:,end:-1:1);
            imageL2(:,:,z) = imread(LFileName2,z)-camOffset;
            tmp = imread(RFileName2,z)-camOffset;
            imageR2(:,:,z) = tmp(:,end:-1:1);
        end
        if t == 0
        % Segment the embryo shell, to compute weights of pixels in fusion.
        level = graythresh(mat2gray(imageL));
        segL   = mat2gray(imageL)>level;
        %
        level  = graythresh(mat2gray(imageR));
        segR   = mat2gray(imageR)>level;
        segR2  = segR;
        segL2  = segL;
        for z = 1 : size(imageL,3)

            res = segL(:,:,z);
            res = bwareaopen(res,1000);
            res = imdilate(res,strel('disk',2));
            res = imfill(res,'holes');
            res = imerode(res,strel('disk',2));
            res = bwareaopen(res,5000);
            %imwrite(res,'segtest.tif','tiff','Compression','None','WriteMode','append');
            segL2(:,:,z) = res;

            res = segR(:,:,z);
            res = bwareaopen(res,1000);
            res = imdilate(res,strel('disk',2));
            res = imfill(res,'holes');
            res = imerode(res,strel('disk',2));
            res = bwareaopen(res,5000);
            %imwrite(res,'segtest.tif','tiff','Compression','None','WriteMode','append');
            segR2(:,:,z) = res;
        end
        disp('Segmentaiton done');
        toc
        % Camera distances and weights; 
        CamDistL = double(cumsum(uint8(segL2),3));
        eta      = 1/90; % weight factor; 
        WeightsL = exp(-eta*CamDistL);
        %WeightsL(:,:,150:end) = 0;
        %WeightsL = WeightsL.*segL2;
        
        CamDistR = double(cumsum(uint8(segR2(:,:,end:-1:1)),3));
        WeightsR = exp(-eta*CamDistR(:,:,end:-1:1));
        %WeightsR(:,:,1:90) = 0;
        %WeightsR = WeightsR.*segR2;
        end
        dynoff  = (median(imageL(segL2>0)))/4+(3*min(imageL(segL2>0)))/4;
        dynoff2 = (median(imageL2(segL2>0)))/4+(3*min(imageL2(segL2>0)))/4;
        
        imageLW = uint16(WeightsL.*double(imageL-dynoff));
        imageRW = uint16(WeightsR.*double(imageR-dynoff));
        
        imageLW2 = uint16(WeightsL.*double(imageL2-dynoff2));
        imageRW2 = uint16(WeightsR.*double(imageR2-dynoff2));
        
        
        disp('Weights done');
        toc

         % rotate weighted right image part ;

         [shiftx,shifty] = xcorr2fft(max(imageRW,[],3),max(imageLW,[],3));
        
        tmp = zeros(size(imageRW,1),size(imageRW,2),...
            size(imageRW,3),class(imageRW));
        
        tmp2 = zeros(size(imageRW2,1),size(imageRW2,2),...
            size(imageRW2,3),class(imageRW2));
        
        if shiftx < 0 
             
            if shifty < 0
                tmp(1:(end+shiftx+1),1:(end+shifty+1),:) = imageRW(-shiftx:end,-shifty:end,:);
                tmp2(1:(end+shiftx+1),1:(end+shifty+1),:) = imageRW2(-shiftx:end,-shifty:end,:);
            else
                tmp(1:(end+shiftx+1),shifty:end,:) = imageRW(-shiftx:end,1:(end-shifty+1),:);
                tmp2(1:(end+shiftx+1),shifty:end,:) = imageRW2(-shiftx:end,1:(end-shifty+1),:);
            end
        else
            
            if shifty < 0
                tmp(max(shiftx,1):end,1:(end+shifty+1),:) = imageRW(1:(end-max(shiftx,1)+1),-shifty:end,:);
                tmp2(max(shiftx,1):end,1:(end+shifty+1),:) = imageRW2(1:(end-max(shiftx,1)+1),-shifty:end,:);
            else
                tmp(max(shiftx,1):end,shifty:end,:) = imageRW(1:(end-max(shiftx,1)+1),1:(end-shifty+1),:);
                tmp2(max(shiftx,1):end,shifty:end,:) = imageRW2(1:(end-max(shiftx,1)+1),1:(end-shifty+1),:);
            end
        end
        
        
%         mX2 = size(imageRW)/2;
%         mX1 = size(imageLW)/2;
%         M = eye(3);
%         
%         T1 = [  1 0 0 -mX2(1)+shifty;
%                 0 1 0 -mX2(2)+shift; 
%                 0 0 1 -mX2(3); 
%                 0 0 0 1 ];
% 
%         RRR = eye(4);
%         RRR(1:3,1:3) = inv(M);
%         % translate back
%         T2 = [  1 0 0 mX1(1);
%                 0 1 0 mX1(2); 
%                 0 0 1 mX1(3); 
%                 0 0 0 1 ];
%         % create a transformation map
%         
%         rotation = diffgeometry.AffineMap(stack.image.domain, stack.image.domain, T2*RRR*T1);
% 
%         rotated = stack.transform(rotation);
%         tmp = rotated.image.apply;
%         disp('rotation done');
%         toc
%         
%         % add up the weigthed images; 
        imout = tmp+imageLW;
        if nch == 2
        imout2 = tmp2+imageLW2;
        end
        for z = 1 : size(imout,3)
            imwrite(imout(:,:,z)',[RFileName(1:end-4),'_f.tif'],'tiff','Compression','none','WriteMode','append');
            if nch == 2
            imwrite(imout2(:,:,z)',[RFileName2(1:end-4),'_f.tif'],'tiff','Compression','none','WriteMode','append');
            end
            %imwrite(imageLW(:,:,z),['fused/',RFileName(1:end-4),'_f_l.tif'],'tiff','Compression','none','WriteMode','append');
            %imwrite(tmp{1}(:,:,z),['fused/',RFileName(1:end-4),'_f_r.tif'],'tiff','Compression','none','WriteMode','append');
        end
        disp('done writeing')
        toc
end