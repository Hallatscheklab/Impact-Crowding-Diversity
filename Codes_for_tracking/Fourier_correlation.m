function [im2_modif]=POC(im1,im2,im2_full)
%fft is the fastest when length of vector is power of 2.
height=size(im1,1);
width=size(im1,2);

fft_im_ref=fft2(im1);
fft_im=fft2(im2);
correl=fftshift((real(ifft2(fft_im_ref.*conj(fft_im)))));
[maxcol,idxline] = max(correl); %max value of each array (row?), and its index in each array (row?)
[max_val,idxcol] = max(maxcol);
dx = idxcol- (width/2+1);
dy = idxline(idxcol)- (height/2+1);
%figure('Name','POC'), imshow(correl,[])
im2_modif = imtranslate(im2_full,[dx, dy]);