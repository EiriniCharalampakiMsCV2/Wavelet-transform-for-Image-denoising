clc 
clear all
close all

% We read the image
im=imread('lena256.bmp');
% returns the double precision value for im.
im=double(im);
[r,c]=size(im);
% Adding noise
im= im+(randn(r,c)*20);
% initialize image, low pass and high pass filters
% low pass filter : smoothing
% high pass filter: edges

low_pass_filter=[0.48296     0.83652     0.22414    -0.12941];

high_pass_filter=[-0.48296 0.83652 -0.22414 -0.12941];

%fliplr is the reverse process

high_pass_filter=fliplr(high_pass_filter);

wavelet_image=zeros(r,c);
temp_image=zeros(r,c);

im_new=imrotate(im,-90);
% we start wavelets, first for rows
for i = 1:r
    % Apply low and high pass filters
    image_after_low_pass = pconv(low_pass_filter,fliplr(im_new(i,:)));
    image_after_high_pass = pconv(high_pass_filter,fliplr(im_new(i,:)));
    
     % Apply downsampling
    image_after_down_sample_low_pass=image_after_low_pass(1:2:length(image_after_low_pass));
    image_after_down_sample_high_pass=image_after_high_pass(1:2:length(image_after_high_pass));

% this is the new image after low+high pass filters and downsampling on
% rows
    temp_image(i,:)=[image_after_down_sample_low_pass image_after_down_sample_high_pass];
end  

% we start downsampling, for columns
for i = 1:c
    % Apply low and high pass filters
    image_after_low_pass = pconv(low_pass_filter,temp_image(:,i)');
    image_after_high_pass = pconv(high_pass_filter,temp_image(:,i)');
    
    % Apply downsampling
    image_after_down_sample_low_pass=image_after_low_pass(1:2:length(image_after_low_pass));
    image_after_down_sample_high_pass=image_after_high_pass(1:2:length(image_after_high_pass));
% this is the new image after low+high pass filters and downsampling on
% rows
    wavelet_image(i,:)=[image_after_down_sample_low_pass image_after_down_sample_high_pass];
end   

 % This is the wavelet image after Level 1
im_new_level1=wavelet_image(1:r/2,1:c/2);

% Now we apply for Level 2 and we will repeat the process, but for 1/4 of
% the image
[r1,c1]=size(im_new_level1);
wavelet_image_level2=zeros(r1,c1);
temp_image_level2=zeros(r1,c1);

im_new_level1=imrotate(im_new_level1,-90);

% We apply wavelets, on rows
for i = 1:r1
    
    image_after_low_pass = pconv(low_pass_filter,fliplr(im_new_level1(i,:)));
    image_after_high_pass = pconv(high_pass_filter,fliplr(im_new_level1(i,:)));
    
    image_after_down_sample_low_pass=image_after_low_pass(1:2:length(image_after_low_pass));
    image_after_down_sample_high_pass=image_after_high_pass(1:2:length(image_after_high_pass));
% this is the new image after low+high pass filters and downsampling on
% rows
    temp_image_level2(i,:)=[image_after_down_sample_low_pass image_after_down_sample_high_pass];
end  

% We apply wavelets, on columns
for i = 1:c1
    
    image_after_low_pass = pconv(low_pass_filter,temp_image_level2(:,i)');
    image_after_high_pass = pconv(high_pass_filter,temp_image_level2(:,i)');
    
    image_after_down_sample_low_pass=image_after_low_pass(1:2:length(image_after_low_pass));
    image_after_down_sample_high_pass=image_after_high_pass(1:2:length(image_after_high_pass));

    % this is the new image after low+high pass filters and downsampling on
% rows
    wavelet_image_level2(i,:)=[image_after_down_sample_low_pass image_after_down_sample_high_pass];
end   

% The wavelet of the Level 2
wavelet_level2=wavelet_image;
wavelet_level2(1:r/2,1:c/2)=wavelet_image_level2;


%% the Reconstruction Process with Soft thresholding
wavelet_level2_new=wavelet_level2;
% Estimate the noise level
hf=[wavelet_level2(129:256,1:128) wavelet_level2(129:256,129:256) wavelet_level2(1:128,129:256)];
sigma=median(abs(hf(:)))/0.6745;
threshold=3*sigma;
% Soft thresholding
wavelet_level2=(sign(wavelet_level2).*(abs(wavelet_level2)-threshold)).*((abs(wavelet_level2)>threshold));

% Apply reconstruction for the wavelet, first we take 1/4 of wavelet Level
% 2 and make reconstruction on it

wavelet_level2_reconstruction=wavelet_level2(1:r/2,1:c/2);

[r2,c2]=size(wavelet_level2_reconstruction);

temp_rec_level1=zeros(r2,c2);
wavelet_level1_rec=zeros(r2,c2);

% We start the process on column
for i =1:c2
    wavelet_level2_low=wavelet_level2_reconstruction(1:(c2/2),i);
    wavelet_level2_high=wavelet_level2_reconstruction((c2/2)+1:c2,i);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
    
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

    % Here we apply the high and low pass filter
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    temp_rec_level1(:,i)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

% We repeat the process for the rows

for i =1:r2
    wavelet_level2_low=temp_rec_level1(i,1:(r2/2));
    wavelet_level2_high=temp_rec_level1(i,(r2/2)+1:r2);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
    
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

     % Here we apply the high and low pass filter
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    wavelet_level1_rec(i,:)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

% Wavelet Level 1 reconstructed
wavelet_level1=wavelet_level2;
wavelet_level1(1:r/2,1:c/2)=wavelet_level1_rec;

% Start the reconstruction process for the original image
wavelet_level2_reconstruction=wavelet_level1;

[r2,c2]=size(wavelet_level2_reconstruction);

temp_rec_level1=zeros(r2,c2);
original_image=zeros(r2,c2);

% the reconstruction starts from the columns

for i =1:c2
    wavelet_level2_low=wavelet_level2_reconstruction(1:(c2/2),i);
    wavelet_level2_high=wavelet_level2_reconstruction((c2/2)+1:c2,i);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

    % Apply low+high pass filters
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    temp_rec_level1(:,i)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

 % We continue the process for the rows
for i =1:r2
    wavelet_level2_low=temp_rec_level1(i,1:(r2/2));
    wavelet_level2_high=temp_rec_level1(i,(r2/2)+1:r2);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

    % Apply low+high pass filters
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    original_image(i,:)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end


%% the Reconstruction Process with Hard thresholding
wavelet_level2_hard=wavelet_level2_new;
% Estimate the noise level
hf=[wavelet_level2_hard(129:256,1:128) wavelet_level2_hard(129:256,129:256) wavelet_level2_hard(1:128,129:256)];
sigma=median(abs(hf(:)))/0.6745;
threshold=3*sigma;
% Hard thresholding
wavelet_level2_hard=(wavelet_level2_hard).*((abs(wavelet_level2_hard)>threshold));

% Apply reconstruction for the wavelet, first we take 1/4 of wavelet Level
% 2 and make reconstruction on it

wavelet_level2_reconstruction=wavelet_level2_hard(1:r/2,1:c/2);

[r2,c2]=size(wavelet_level2_reconstruction);

temp_rec_level1=zeros(r2,c2);
wavelet_level1_rec=zeros(r2,c2);

% We start the process on column
for i =1:c2
    wavelet_leavel2_low=wavelet_level2_reconstruction(1:(c2/2),i);
    wavelet_level2_high=wavelet_level2_reconstruction((c2/2)+1:c2,i);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
    
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

    % Here we apply the high and low pass filter
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    temp_rec_level1(:,i)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

% We repeat the process for the rows

for i =1:r2
    wavelet_level2_low=temp_rec_level1(i,1:(r2/2));
    wavelet_level2_high=temp_rec_level1(i,(r2/2)+1:r2);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
    
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

     % Here we apply the high and low pass filter
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    wavelet_level1_rec(i,:)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

% Wavelet Level 1 reconstructed
wavelet_level1_hard=wavelet_level2_hard;
wavelet_level1_hard(1:r/2,1:c/2)=wavelet_level1_rec;

% Start the reconstruction process for the original image
wavelet_level2_reconstruction=wavelet_level1_hard;

[r2,c2]=size(wavelet_level2_reconstruction);

temp_rec_level1=zeros(r2,c2);
original_image_hard=zeros(r2,c2);

% the reconstruction starts from the columns

for i =1:c2
    wavelet_level2_low=wavelet_level2_reconstruction(1:(c2/2),i);
    wavelet_level2_high=wavelet_level2_reconstruction((c2/2)+1:c2,i);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

    % Apply low+high pass filters
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    temp_rec_level1(:,i)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

 % We continue the process for the rows
for i =1:r2
    wavelet_level2_low=temp_rec_level1(i,1:(r2/2));
    wavelet_level2_high=temp_rec_level1(i,(r2/2)+1:r2);

    wavelet_level2_up_low=zeros(1,2*length(wavelet_level2_low));
    wavelet_level2_up_high=zeros(1,2*length(wavelet_level2_high));
% upsampling
    wavelet_level2_up_low(1:2:length(wavelet_level2_up_low))=wavelet_level2_low;
    wavelet_level2_up_high(1:2:length(wavelet_level2_up_high))=wavelet_level2_high;

    % Apply low+high pass filters
    wavelet_level2_up_low=pconv(low_pass_filter,fliplr(wavelet_level2_up_low));
    wavelet_level2_up_high=pconv(high_pass_filter,fliplr(wavelet_level2_up_high));
    
    original_image_hard(i,:)=fliplr(wavelet_level2_up_low+wavelet_level2_up_high);
end

%% display the results

figure(1)
subplot(3,3,1)
imshow(im,[])
title('The original image')
subplot(3,3,2)
imshow(wavelet_image,[])
title('Level 1 wavelet')
subplot(3,3,3)
imshow(wavelet_level2,[])
title('Level 2 wavelet')
subplot(3,3,6)
imshow(wavelet_level2_hard,[])
title('Level 2 wavelet with soft threshold')
subplot(3,3,5)
imshow(wavelet_level1,[])
title('reconstructed Level 1 wavelet with soft threshold')
subplot(3,3,4)
imshow(original_image,[])
title('reconstructed original image with soft threshold')
subplot(3,3,9)
imshow(wavelet_level2_hard,[])
title('Level 2 wavelet with hard threshold')
subplot(3,3,8)
imshow(wavelet_level1_hard,[])
title('reconstructed Level 1 wavelet with hard threshold')
subplot(3,3,7)
imshow(original_image_hard,[])
title('reconstructed original image with hard threshold')
