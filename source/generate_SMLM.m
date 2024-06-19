%% This is code is for generating image-stack of SMLM raw data

Nlow = 6000;  % lower bound of total photons
Nup = 6000; % upper bound of total photons
bg = 100; % photons of background noise 
Nmid = (Nlow+Nup)/2;
Nr = Nup - Nmid;
total_round = 8000;  % total frames

%% 0. Determine parameters and read relevant data
for i = 1:size(real_PSF,3) % Normalize
    tmp = real_PSF(:,:,i);
    tmp = tmp/sum(tmp(:));
    real_PSF(:,:,i) = tmp;
end
coeff = Spline3D_interp(real_PSF);


off = 3;
[x, y, z , ~] = size(coeff);
xsize = x+1;
zsize = z;  
zz = (zsize)/2;
ROI = xsize - 2*off;
r_ROI = (ROI-1)/2;

% read the real location of simulated biological structures
data = readmatrix('csv data\EPFL.csv');
x_real = data(:,3);  
y_real = data(:,4);
z_real = data(:,5);


num = 1;
range1 = zeros(3,3);
min_list = [min(x_real),min(y_real),min(z_real)];
max_list = [max(x_real),max(y_real),max(z_real)];
range1(1,:) = min_list;
range1(2,:) = max_list;
range1(3,:) = max_list - min_list;

x_real_xpixel = x_real/x_pixelsize;
y_real_xpixel = y_real/x_pixelsize;
z_real_zpixel = z_real/z_pixelsize;

% Perform mapping to generate a rough image
show_pixel = 10;  % The pixel size of the final displayed super-resolution image
img_full = show_full(show_pixel,x_real, y_real);
figure;imshow(img_full.^0.3,[]);title('img full')

x_real_range = [min(x_real) max(x_real)]
z_real_range = [min(z_real) max(z_real)]

x_real_pixel_range = [min(x_real_xpixel) max(x_real_xpixel)]
z_real_pixel_range = [min(z_real_zpixel) max(z_real_zpixel)]

%% 1. Start simulation
x_final1 = x_real_xpixel;
y_final1 = y_real_xpixel;
z_final1 = z_real_zpixel;


% region selection
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
xhalf = (7500)/x_pixelsize;  % unit pixel
% xhalf = (2000)/x_pixelsize;  % unit pixel
yhalf = xhalf;

% EPFL.csv
xcenter = (floor(size(img_full,1)/2)+3)*show_pixel/x_pixelsize; % xcenter和ycenter适度修改  
ycenter = (floor(size(img_full,2)/2)+3)*show_pixel/x_pixelsize;


% % NPC0905(_2).csv 
% xcenter = (floor(size(img_full,1)/2)+3)*show_pixel/x_pixelsize; % xcenter和ycenter适度修改  
% ycenter = (floor(size(img_full,2)/2)+3)*show_pixel/x_pixelsize;

% data1.csv 
% xcenter = (1346)*show_pixel/x_pixelsize; % xcenter和ycenter适度修改 
% xcenter = (500)*show_pixel/x_pixelsize;
% ycenter = (500)*show_pixel/x_pixelsize;


x_low = xcenter - xhalf + min(x_final1);
x_up = xcenter + xhalf + min(x_final1);
y_low = ycenter - yhalf + min(y_final1);
y_up = ycenter + yhalf + min(y_final1);
mask = x_final1>x_low & x_final1<x_up & y_final1>y_low & y_final1<y_up;
x_final1 = x_final1(mask);
y_final1 = y_final1(mask);
z_final1 = z_final1(mask);

%==============================================================================================================================================================================
x_final1 = x_final1-min(x_final1);
y_final1 = y_final1-min(y_final1);
%==============================================================================================================================================================================

x_final = repmat( x_final1,num,1);
y_final = repmat( y_final1,num,1);
z_final = repmat( z_final1,num,1);
numall = size(x_final1,1); 
range_z = [min(z_final1)*z_pixelsize, max(z_final1)*z_pixelsize]  
locs = numall*num  

img_fullxy = show_full(show_pixel,x_final1*x_pixelsize, y_final1*x_pixelsize);
img_fullxz = show_full(show_pixel,x_final1*x_pixelsize, z_final1*z_pixelsize);
img_fullyz= show_full(show_pixel,y_final1*x_pixelsize, z_final1*z_pixelsize);


figure;imshow(img_fullxy.^0.3,[]);title('img ROI(xy)') %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
figure;imshow(img_fullxz.^0.3,[]);title('img ROI(xz)') %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
figure;imshow(img_fullyz.^0.3,[]);title('img ROI(yz)') %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





close all
%% 2. Specify the parameters for each subgraph
% 50*50 μm 140 points
dis = 20;
rangex = (max(x_final1) - min(x_final1));
rangey = (max(y_final1) - min(y_final1));
max_range = max([rangex rangey]);
if max_range < 2*xhalf
    final_range = max_range;
    imagesize = floor(max_range)+1+2*(r_ROI+dis);
else
    final_range = 2*xhalf;
    imagesize = floor(2*xhalf)+1+2*(r_ROI+dis); 
end
density = floor(140 * (max_range/500)^2)+2;   % The number of blinking points in each subgraph

N = rand(density,total_round);
N = (2*N-1)*Nr + Nmid;

image_stack = zeros(imagesize, imagesize, total_round);
randmatrix = randi([1 locs],density,total_round);

x_rough = round(x_final);   
y_rough = round(y_final);   
x_loc = x_final - x_rough;   
y_loc = y_final - y_rough;   
z_loc = z_final;

tic
parfor index1 = 1:total_round
    tmp_image = zeros(imagesize,imagesize);
    tmp_xloc = x_loc(randmatrix(:,index1),:);
    tmp_xrough = x_rough(randmatrix(:,index1),:)+(r_ROI+dis);
    tmp_yloc = y_loc(randmatrix(:,index1),:);
    tmp_yrough = y_rough(randmatrix(:,index1),:)+(r_ROI+dis);
    tmp_zloc = z_loc(randmatrix(:,index1),:);
    tmp_Nlist = N(:,index1);

    for index = 1:density
        xc = tmp_xloc(index);
        yc = tmp_yloc(index);
        zc = zz + tmp_zloc(index);

        delta_x = -1*xc;
        xstart = floor(delta_x);
        delta_x = delta_x - xstart;

        delta_y = -1*yc;
        ystart = floor(delta_y);
        delta_y = delta_y - ystart;

        delta_z = zc - floor(zc);
        [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2((delta_x),(delta_y),(delta_z));

        z_index = floor(zc);

        if z_index <=0
            z_index = 1;
        end
        if z_index >z
            z_index = z;
        end

        test1 = zeros(ROI, ROI);
        for ii = 0:ROI-1
            for jj = 0:ROI-1
                temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,z,delta_f,coeff);
                test1(ii+1,jj+1)=temp;
            end
        end
        test1 = test1/sum(test1(:));
        tmp_image(tmp_xrough(index)-r_ROI:tmp_xrough(index)+r_ROI, tmp_yrough(index)-r_ROI:tmp_yrough(index)+r_ROI)=...
            tmp_image(tmp_xrough(index)-r_ROI:tmp_xrough(index)+r_ROI, tmp_yrough(index)-r_ROI:tmp_yrough(index)+r_ROI)+tmp_Nlist(index)*test1;
        
    end

    
    tmp_image = uint16(tmp_image) + bg;   % add background noise
    tmp_image = imnoise(tmp_image, 'poisson');   % add shot noise (poisson noise)
    image_stack( :, :, index1) = double(tmp_image);
%     figure;imshow(tmp_image,[])

end
toc

ims = image_stack(:,:,1:6000);

imageslicer(single(permute(ims(:,:,1:100),[2 1 3])));
title('SMLM raw data(Fisrt 100 frames)')














