function imagebgc = TopHatFFT(image,aspect,BallRad,IBd)

    if nargin == 1
        
        BallRad = 100; % radius of the ball for background estimation
        IBd = [0, 2^16-1]; % boundaries of intensities
        aspect = 1.25/.26;% aspect ratio in image;
    end
    
    % construct the structural element
    se = zeros([2*BallRad+1 2*BallRad+1 2*BallRad+1]);
    [X,Y,Z] = meshgrid(-BallRad:BallRad, -BallRad:BallRad, -BallRad:BallRad);
    se(X.^2 + Y.^2 + aspect.^2*Z.^2 <= BallRad^2 + 1) = 1;
    clear X Y Z 
    se = se/sum(se(:));
    seBig = zeros(size(image));
    seBig(round(end/2-BallRad):round(end/2+BallRad),round(end/2-BallRad):round(end/2+BallRad),round(end/2-BallRad):round(end/2+BallRad)) = se;
    clear se
    
    % scale image intensities to [0 1]
    img = mat2gray(image,IBd);

    % erosion is a dilation on the complementary image; 
    imc = 1-img;

    % convolution theorem.
    imageft = fftn(img);
    seBigft = fftn(seBig);
    conft = imageft.*seBigft;
    
    % get eroded image by using the complementary image of the inverse fft.
    erode = 1-fftshift(ifftn(conft));
    
    % now dilate this, to finish the tophat
    erodeft = fftn(erode);
    conft = erodeft.*seBigft;
    
    con = 1-fftshift(ifftn(conft));
    
    imagebgc = img-con;
    imagebgc(imagebgc<0) = 0;
    
    imagebgc = uint16(imagebgc*IBd(2));
end

