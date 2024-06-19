% Script for calculating image shift
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
% 
function [shift1, shift2, similarity]=registration_in_each_channel(im1,im2)
% [shift1, shift2, tmpval1] = registration_in_each_channel(ref_plane1(:,:,jj), img_plane1(:,:,ii));
% im1是reference image(大小为32*32,double，是某一个深度的reference)；
% im2是img_plane1(大小为32*32 pixels single)，表示的是


output = dftregistration(fft2(im1),fft2(im2),10);   % 这个是subpixel registration的方法

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
%
% Outputs
% output =  [error, diffphase, net_row_shift, net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.



shift1 = output(3);  % net_row_shift
shift2 = output(4);  % net_col_shift

im2_shift = FourierShift2D(im2, [shift1 shift2]);     % 这是针对每一个ROI的移到和reference相同的位置，在仔细看看

[similarity,cc] = cc2(im1(3:end-2,3:end-2),im2_shift(3:end-2,3:end-2));