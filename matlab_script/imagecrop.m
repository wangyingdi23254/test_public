%%
clc;
clear all;
path = pwd;
%image file direction
dirname1 = [path,'\T'];
dirname2 = [path,'\rou'];
mkdir(dirname1); 
mkdir(dirname2);
imfiles = dir([path '\*.tif']);% original images input
imnum = length(imfiles);
for i=1:imnum
    imname = imfiles(i).name;
    im = imread([path '\' imname]);
    %[x y]=size(im);
    im1 = imcrop(im,[0,0,1072,630]); % crop image by rectangle [x1,y1,x2,y2]
    imwrite(uint8(im1),[dirname1, '\', num2str(i,'%03d'),'.tif'],'tif')
    im2 = imcrop(im,[0,630,2072,1260]);
    imwrite(uint8(im2),[dirname2, '\', num2str(i,'%03d'),'.tif'],'tif')
end

