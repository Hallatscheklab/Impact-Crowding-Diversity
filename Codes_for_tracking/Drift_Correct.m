currentpath=pwd;
filename = 'XXX.tif'; %Enter the file name.
path = fullfile(currentpath,filename);
info = imfinfo(path);
MaxIndex = length(info);

%size of image should be power of 2, for rapid FFT.

% set reference pic and select ROI 
im1_full=imread(path,'index',1);
figure('Name','ROI_select'), imshow(im1_full);
ROI=getrect();
xmin  =ROI(1);
ymin  =ROI(2);
width =ROI(3);
height=ROI(4);
xcenter = xmin + width/2.0;
ycenter = ymin + height/2.0;
crop_width=2^(nextpow2(double(width)))-1;
crop_height=2^(nextpow2(double(height)))-1;
im1_ROI=im1_full(floor(ycenter-crop_height/2.0):floor(ycenter+crop_height/2.0),floor(xcenter-crop_width/2.0):floor(xcenter+crop_width/2.0));
figure('Name','selected_ROI'), imshow(im1_ROI);

%make directory and save first pic
mkdir(fullfile(currentpath,'Drift-corrected'));
imwrite(uint8(im1_full),fullfile(currentpath,'Drift-corrected','Drift_F.tif'));

%Drift correction
for i = 2 : MaxIndex
    im2_full = imread(path,'index',i);
    im2_ROI=im2_full(floor(ycenter-crop_height/2.0):floor(ycenter+crop_height/2.0),floor(xcenter-crop_width/2.0):floor(xcenter+crop_width/2.0));
    
    immodif2 = Fourier_correlation(im1_ROI,im2_ROI,im2_full);
    imwrite(uint8(immodif2),fullfile(currentpath,'Drift-corrected','Drift_F.tif'),'WriteMode','append');
end
