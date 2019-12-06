nameDummy = 'TimeCol_2_MMStack';
name = [nameDummy,'.ome.tif'];
image = imread(name);
crop_X1 = 1:size(image,2);
crop_X2 = crop_X1;
crop_Y1 = 1:size(image,1);
crop_Y2 = crop_Y1;

RealStackSize = 451;
TimeCount = 0;
counter = 1;

col = 1;
ncol = 2;

image = zeros(length(crop_Y1),length(crop_X1),RealStackSize*2,class(image));
%image1 = zerzs(size(image,1),length(crop_X1),RealStackSize,class(image));
%image2 = zeros(size(image,1),length(crop_X2),RealStackSize,class(image));
%%
for time =  0:24
    %name = ['Time_1_MMStack_',sprintf('%01d',time),'.ome.tif'];
    
 %   for z = 1 : RealStackSize
    if time == 0
        
    else   
        name = [nameDummy,'_',num2str(time),'.ome.tif'];
    end
    StackSize = length(imfinfo(name));
    
  
    for k = 1 : StackSize
        
        %image(:,:,counter) = imread(name,k);
        temp  = imread(name,k);
        image(:,:,counter) = temp(crop_Y1,crop_X1);
        %image2(:,:,counter) = temp(:,crop_X2);
        
        counter = counter+1;
        
        if counter > RealStackSize*2
         
            
           counter = 1;

           name1  =  ['Time_',sprintf('%06d',TimeCount),'.tif'];
           for j = 1 : RealStackSize/1
               imwrite(image(:,:,2*j),[name1(1:end-4),'_left_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
               imwrite(image(:,:,2*(j-1)+1),[name1(1:end-4),'_right_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
           end
           imwrite(max(image(:,:,2:2:end),[],3),[name1(1:end-4),'_max_left_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');
           imwrite(max(image(:,:,1:2:end),[],3),[name1(1:end-4),'_max_right_crop_c',num2str(col),'.tif'],'Compression','none','WriteMode','Append');

                          
           col = col+1;
           if col > ncol
               TimeCount = TimeCount+1;
               col = 1;
           end
           
        end
        
        
        
        
    end
    
end

