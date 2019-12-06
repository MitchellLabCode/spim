crop_X1 = 1:1204;
crop_X2 = crop_X1+1204;

RealStackSize = 411;
TimeCount = 80;
counter = 1;
%%
nameDummy = 'MyosinOnly00034_P';
name = [nameDummy,num2str(0),'.TIF'];
image = imread(name);
%image = zeros(size(image,1),size(image,2),RealStackSize,class(image));
image1 = zeros(size(image,1),length(crop_X1),RealStackSize,class(image));
image2 = zeros(size(image,1),length(crop_X2),RealStackSize,class(image));
%%
for time = 86:89
    %name = ['Time_1_MMStack_',sprintf('%01d',time),'.ome.tif'];
    
 %   for z = 1 : RealStackSize
    
    name = [nameDummy,num2str(time),'.TIF'];
    StackSize = length(imfinfo(name));
    
  
    for k = 1 : StackSize
        
        
        if counter > RealStackSize
         
            
           counter = 1;

           name1  =  ['Time_',sprintf('%06d',TimeCount),'.tif'];
           for j = 1 : RealStackSize/1
               imwrite(image1(:,:,j),[name1(1:end-4),'_left_crop.tif'],'Compression','none','WriteMode','Append');
               imwrite(image2(:,:,j),[name1(1:end-4),'_right_crop.tif'],'Compression','none','WriteMode','Append');
           end
           
           TimeCount = TimeCount+1;
        end
        
        %image(:,:,counter) = imread(name,k);
        temp  = imread(name,k);
        image1(:,:,counter) = temp(:,crop_X1);
        image2(:,:,counter) = temp(:,crop_X2);
        
        counter = counter+1;
    end
    
end

