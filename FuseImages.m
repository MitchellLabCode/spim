function [imout,shiftx,shifty] = FuseImages(image1,image2,weight1,weight2)

    % weight is an optional input, should contain to weight images,
    % specifying the contribution of image1 and image2 to the final fused
    % image;
    
    
    [shiftx,shifty] = xcorr2fft(max(image1,[],3),max(image2,[],3));
    
    image1 = uint16(weight1.*double(image1));
    image2 = uint16(weight2.*double(image2));
    
    tmp = zeros(size(image1,1),size(image1,2),...
    size(image1,3),class(image1));


    if shiftx < 0 

        if shifty < 0
            tmp(1:(end+shiftx+1),1:(end+shifty+1),:) = image1(-shiftx:end,-shifty:end,:);
        else
            tmp(1:(end+shiftx+1),max(shifty,1):end,:) = image1(-shiftx:end,1:(end-max(shifty,1)+1),:);
        end
    else

        if shifty < 0
            tmp(max(shiftx,1):end,1:(end+shifty+1),:) = image1(1:(end-max(shiftx,1)+1),-shifty:end,:);
        else
            tmp(max(shiftx,1):end,max(shifty,1):end,:) = image1(1:(end-max(shiftx,1)+1),1:(end-max(shifty,1)+1),:);
        end
    end

    imout = tmp+image2;
end
