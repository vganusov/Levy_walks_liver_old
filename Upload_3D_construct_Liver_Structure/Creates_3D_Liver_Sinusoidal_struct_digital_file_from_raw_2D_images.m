% Creates the binary data file of the 3D structure  

clear all; close all
clf
k=1;
for i=1:21;%21;%20;
for j=1:2;
['17_1207_evans_blue_10_percent_liver_2_0002_',num2str(i)]
imgT(:,:,:) = imread(['17_1207_evans_blue_10_percent_liver_2_0002_',num2str(i),'.tif']);
img(:,:,:)=imgT(:,:,1:3) ; 
figure()
imshow(img);
img2 = rgb2gray(img(:,:,1:3));%round(img(1:10,1:10)./100));
% remove some noise
img_gauss = imfilter(img2,fspecial('gaussian',10,2));%10,2
 %figure();
 %imshow(img_gauss); 
% remove all the lines from the image
img_filtered = medfilt2(img_gauss,[31,31],'symmetric');
% figure();
% imshow(img_filtered);
% subtract image without lines from image with lines
img_subtracted = uint8(double(img_gauss)-double(img_filtered)+127);
% figure();
% imshow(img_subtracted);
% threshold
img_thresh = img_subtracted>130;%130
% figure();
% imshow(img_thresh);
% use open to remove some noise
%figure();
img_opened = imopen(img_thresh,ones(3));
imshow(img_opened);
IMG_OPENED_3D_Evans_old(:,:,k)=img_opened;
k=k+1;
%figure()
BW3 = bwmorph(img_opened,'skel',Inf);
%imshow(BW3)
%BW3_3D(:,:,i)=BW3;
end
end
save IMG_OPENED_3D
