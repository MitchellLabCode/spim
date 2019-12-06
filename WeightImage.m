function [Weights,seg2] = WeightImage(image,camFlag,eta,disk,amin,amin2)

    if nargin < 3
        eta = 1/35;
        disk = 2;
        amin = 1000;
        amin2 = 5000;
    end
    image = image;
    % Segment the embryo shell, to compute weights of pixels in fusion.
    level = graythresh(mat2gray(image));
    seg   = mat2gray(image)>level;
    seg2  = seg;
    for z = 1 : size(image,3)

        res = seg(:,:,z);
        res = bwareaopen(res,amin);
        res = imdilate(res,strel('disk',disk));
        res = imfill(res,'holes');
        res = imerode(res,strel('disk',disk));
        res = bwareaopen(res,amin2);
        seg2(:,:,z) = res;
    end
    disp('Segmentation done');
    
    %eta      = 1/90; % weight factor; 
    if camFlag == 1
        % Camera distances and weights; 
        CamDist = double(cumsum(uint8(seg2(:,:,end:-1:1)),3));
        Weights = exp(-eta*CamDist(:,:,end:-1:1));
    else
        % Camera distances and weights; 
        CamDist = double(cumsum(uint8(seg2),3));
        Weights = exp(-eta*CamDist);
    end