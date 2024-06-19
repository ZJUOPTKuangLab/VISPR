close all;
clear
clc

% We choose different zernike number to experiment. In most cases, we
% choose 21, 32, 45, 60 items for comparison. We found that the more items
% there are, the better the fitting effect is, but is more sensitive to noise. 
% Therefore, in general experiments, we choose 21 items as the default value.

global zernikenum;
zernikenum = 21;


%% 1. Generate ground truth vector PSF
addpath('source\')
addpath('utilities\')
addpath('shared\')

parameters_real.NA = 1.45;    
parameters_real.refmed = 1.518;   % sample fluid refraction
parameters_real.refcov = 1.518;   % cover glass refraction
parameters_real.refimm = 1.518;   % immersion oil refraction
% parameters_real.NA = min(parameters_real.NA, parameters_real.refmed);
parameters_real.lambda = 590;
parameters_real.pixelSizeX = 89;       % nm, pixel size of the image
parameters_real.pixelSizeY = 89;       % nm, pixel size of the image

parameters_real.zemit = -1000:10:1000;     %nm z positions of each PSF
Nz = numel(parameters_real.zemit);
parameters_real.sizeZ = Nz;              % number of discrete z positions
parameters_real.xemit = 0;                %nm x positions of each PSF
parameters_real.yemit = 0;                %nm y positions of each PSF
parameters_real.objStage = parameters_real.zemit*0;   %nm objective positions of each PSF
parameters_real.deltaz = parameters_real.zemit(1,end)-parameters_real.zemit(1,end-1);
parameters_real.objStage0 = 0;              % nm, distance from coverslip to nominal focal plane, when RI match,this is useless
parameters_real.zemit0 = -1*parameters_real.refmed/parameters_real.refimm*(parameters_real.objStage0); % reference emitter z position, nm, distance of molecule to coverslip
parameters_real.ztype = 'emitter';

Npixels = 31;                           % number of pixels of the psf image
parameters_real.sizeX = Npixels;             % number of pixels of the psf image
parameters_real.sizeY = Npixels;             % number of pixels of the psf image
parameters_real.Npupil = 64;                 % sampling at the pupil plane

% total zernike aberration orders, n,m,amplitude(nm)
parameters_real.aberrations = [2		2	60
                                2		-2	-20
                                3		1	30
                                3		-1	0
                                4		0	20
                                3		3	31
                                3		-3	0
                                4		2	0
                                4		-2	0
                                5		1	0
                                5		-1	0
                                6		0	15
                                4		4	0
                                4		-4	0
                                5		3	0
                                5		-3	0
                                6		2	0
                                6		-2	0
                                7		1	0
                                7		-1	0
                                8		0	0];
%                                 5	    5	0
%                                 5	    -5	0
%                                 6	    4	0
%                                 6	    -4	0
%                                 7	    3	0
%                                 7	    -3	0
%                                 8	    2	0
%                                 8	    -2	0
%                                 9	    1	0
%                                 9	    -1	0
%                                 10	    0	0];
%                                 6	    6	0
%                                 6	    -6	0
%                                 7	    5	0
%                                 7	    -5	0
%                                 8	    4	0
%                                 8	    -4	0
%                                 9	    3	0
%                                 9	    -3	0
%                                 10	    2	0
%                                 10	    -2	0
%                                 11	    1	0
%                                 11	    -1	0
%                                 12	    0	0
%                                 7	    7	0
%                                 7	    -7	0
%                                 8	    6	0
%                                 8	    -6	0
%                                 9	    5	0
%                                 9	    -5	0
%                                 10	    4	0
%                                 10	    -4	0
%                                 11	    3	0
%                                 11	    -3	0
%                                 12	    2	0
%                                 12	    -2	0
%                                 13	    1	0
%                                 13	    -1	0
%                                 14	    0	0];

Nphotons = 2000 + 0*ones(1,Nz);    % photons for each PSF
bg = 0+ 0*ones(1,Nz);          % background for each PSF
parameters_real.Nphotons = Nphotons;
parameters_real.bg = bg;

% Bead parameters for convolution with PSF and derivatives, beaddiameter in nm
parameters_real.bead = false;
parameters_real.beaddiameter = 100;
% check on meaningfullness of bead convolution
if parameters_real.beaddiameter<=parameters_real.pixelSizeX
  parameters_real.bead = false;
end
% generate real vectorial PSF
[real_PSF] = gen_vectorial_psf(parameters_real);
imageslicer(single(permute(real_PSF(:,:,1:10:end),[2 1 3])));
title('real PSF(interval 100nm)')

% Draw the corresponding Pupil and Zernike coefficients for PSF
pupil_real= Zernike_construct_pupil(parameters_real);
figure;imagesc(pupil_real);title('Zernike pupil(GT)');drawnow;
f1 = figure;
axMode =axes(f1);
bar(parameters_real.aberrations(:,3),0.3);
legend(axMode,'Real Zernike Coefficient')
orders = parameters_real.aberrations(:,1:2);
for k=zernikenum:-1:1
    axn{k}=[num2str(orders(k,1)) ',' num2str(orders(k,2))];
end
axMode.XTick=1:length(axn);
axMode.XTickLabel=axn;
ylabel('Zernike rms value(nm)') ,xlabel('Zernike mode')
%% 2. Generate simulation single-molecule dataset
x_pixelsize = parameters_real.pixelSizeX;
z_pixelsize = parameters_real.zemit(2) - parameters_real.zemit(1);
run generate_SMLM.m

%% 3. segmentation
setup.is_imgsz = 1 ;
setup.is_sCMOS = 0;
setup.offset = 0;
setup.gain = 1;

boxsz = 31;
thresh_dist = 20;
thresh = 1.5; % threshold needs to be carefully selected, otherwise the segmentation effect will not be ideal
setup.use_default_thresh = 0;
[subregion_ch1, seg_display] = crop_subregion_VISPR(ims, boxsz, thresh_dist, thresh, setup);

% boxsz = 25;
% thresh_dist = 20;
% thresh_low = 50;
% thresh_high = 55;
% thresh = [thresh_low thresh_high]; % threshold needs to be carefully selected, otherwise the segmentation effect will not be ideal
% [subregion_ch1, seg_display] = crop_subregion_ast(ims,boxsz,thresh,thresh_dist,setup);

num_subregion = size(subregion_ch1,3);
disp(['Finish segmentation! There are ' num2str(num_subregion) ' subregion images']);
total_point = size(subregion_ch1,3);
num_display = 3;  % frame index 
if num_display < 1 || num_display > size(seg_display.ims_ch1,3)
    msgbox('Please input the correct number!');
end
disp('Show segmentation results');
raw = seg_display.ims_ch1(:,:,num_display)/max(max(seg_display.ims_ch1(:,:,num_display)));
index_rec = find(seg_display.allcds_mask(:,3) == num_display-1);
rec_vector = cat(2,seg_display.t1(index_rec),seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
figure; imshow(img_select);
axis tight
title('Segmentation results');
subregion_ch1 = subregion_ch1(:,:,randi(total_point,1,min([total_point,5000])));
% Image preprocessing: subtract the background
% mmed = quantile(subregion_ch1(:),0.01);
% imbg = subregion_ch1(subregion_ch1<mmed);
% BG =  mean(imbg(:));
% subregion_ch1 = subregion_ch1 - BG;
% subregion_ch1(subregion_ch1<0) = 1e-6;
% subregion_ch1 = abs(subregion_ch1 - BG);
for i = 1: size(subregion_ch1,3)
    tmp = subregion_ch1(:,:,i);
    BG = min(tmp(:));
    tmp = tmp - BG;
    subregion_ch1(:,:,i) = tmp;
end
imageslicer(subregion_ch1);
% Obtain preprocessed background and photons
mmed = quantile(subregion_ch1(:),0.3);
imbg = subregion_ch1(subregion_ch1<mmed);
bg =  mean(imbg(:));
Nph = mean(sum(sum(subregion_ch1-bg)));
Npixels = size(subregion_ch1,1);

%% VISPR to get accurate 3D insitu PSF model
zemit = -800:50:800;
iteration = 8;
empupil.min_similarity = 0.8;
empupil.bin_lowerBound = 50;
Zernike = parameters_real.aberrations;
standard = 'Wyant';
if parameters_real.bead == true
    run VISPR_beads.m
else
    run VISPR.m
end
% Generating corresponding 3D insitu PSF model
parameters_insitu = set_parameters(parameters_real.sizeX, numel(parameters_real.zemit), parameters_real.zemit,standard);
parameters_insitu.NA = parameters_real.NA;
parameters_insitu.refmed = parameters_real.refmed;
parameters_insitu.refcov = parameters_real.refcov;
parameters_insitu.refimm = parameters_real.refimm;
parameters_insitu.lambda = parameters_real.lambda;
parameters_insitu.pixelSizeX = parameters_real.pixelSizeX;
parameters_insitu.pixelSizeY = parameters_real.pixelSizeY;
parameters_insitu.aberrations(:,3) = P(1:zernikenum);
insitu_PSF = gen_vectorial_psf(parameters_insitu);

%% Compare real PSF with insitu PSF
imageslicer(single(permute(real_PSF(:,:,1:10:end),[2 1 3])));
title('real PSF(interval 100nm)')

imageslicer(single(permute(insitu_PSF(:,:,1:10:end),[2 1 3])));
title('insitu PSF(interval 100nm)')

%% Calculatie the Localization Precision and Accuracy
for i = 1:size(real_PSF,3)
    tmp1 = real_PSF(:,:,i);
    real_PSF(:,:,i)=tmp1/sum(tmp1(:));

    tmp2 = insitu_PSF(:,:,i);
    insitu_PSF(:,:,i)=tmp2/sum(tmp2(:));
end
coeff_real = Spline3D_interp(real_PSF);
coeff_VISPR = Spline3D_interp(insitu_PSF);

run calculate_preformance.m
%% 6. Reconstrction of simulated  microtubule
boxsz1 = 19;
thresh_dist1 = 14;
thresh = 1.5;
[subregion1, seg_display1] = crop_subregion_VISPR(image_stack, boxsz1,thresh_dist1, thresh, setup);

num_subregion1 = size(subregion1,3);
disp(['Finish segmentation! There are ' num2str(num_subregion1) ' subregion images']);
x_list = seg_display1.allcds_mask(:,2);
y_list = seg_display1.allcds_mask(:,1);
x_list = x_list+5; 
y_list = y_list+5; % Compenstate for the segmentation
frame = seg_display1.allcds_mask(:,3);
r_boxsz = (boxsz1-1)/2;
x_start = single(x_list- r_boxsz);
y_start = single(y_list- r_boxsz);
total_candidates = num_subregion1;

batchsize = 500000; 
if mod(num_subregion1,batchsize)==0
    it_all = num_subregion1/batchsize;
else
    it_all = floor(num_subregion1/batchsize)+1;
end
sCMOSvarmap = 0;
subregion1 = single(subregion1);
Pcspline_real = zeros(num_subregion1, 6);
Pcspline_XX = zeros(num_subregion1, 6);
Pcspline_HF = zeros(num_subregion1, 6);
uncertainty_real = zeros(num_subregion1, 5);
loc_ll_real = zeros(num_subregion1, 1);

% 6.2 real PSF Reconstruction =========================================================================================================================================================================
for seq = 1:it_all
    if seq == it_all
        range = batchsize*(seq-1)+1:total_candidates;
    else
        range = batchsize*(seq-1)+1:batchsize*seq;
    end
    tic
    [t_Pcspline,CRLB,LogL]=mleFit_LM(subregion1(:, :, range), 6, 50,single(coeff_real),sCMOSvarmap,1);  % 
    toc
%     Pcspline(range, :) =  [t_Pcspline,CRLB,LogL];
    Pcspline_real(range, :) =  t_Pcspline; % x,y,N,bg,z
    uncertainty_real(range,:) = CRLB; % x,y,N,bg,z
    loc_ll_real(range,:) = LogL;
end

% remove low confidence spots
result_real = zeros(num_subregion1,5); % N, x, y, z, bg
result_real(:,1) = Pcspline_real(:,3);
result_real(:,2) = Pcspline_real(:,1);  % x
result_real(:,3) = Pcspline_real(:,2);  % y
result_real(:,4) = Pcspline_real(:,5);  % z
result_real(:,5) = Pcspline_real(:,4);
result_real = single(result_real);

result_uncertainty_real = zeros(num_subregion1,5);
result_uncertainty_real(:,1) = uncertainty_real(:,3); % N, x, y, z, bg
result_uncertainty_real(:,2) = uncertainty_real(:,1);
result_uncertainty_real(:,3) = uncertainty_real(:,2);
result_uncertainty_real(:,4) = uncertainty_real(:,5);
result_uncertainty_real(:,5) = uncertainty_real(:,4);

loc_crlbx = result_uncertainty_real(:,2);
loc_crlby = result_uncertainty_real(:,3);
loc_crlbz = result_uncertainty_real(:,4);

loc_photons = result_real(:,1); 
loc_x = x_start + result_real(:,2); % pixel
loc_y = y_start + result_real(:,3); % pixel
loc_z = result_real(:,4)- size(coeff,3)/2; % pixel

llthreshold = -600;
min_photon = 1000;
loc_uncer_max = (1/x_pixelsize);
loc_uncer_zmax = (1);  
zmask_high = 800/z_pixelsize;
zmask_low = -800/z_pixelsize;
llmask = loc_ll_real > llthreshold;
intmask = loc_photons < min_photon;
uncermask = (loc_crlbx > loc_uncer_max) | (loc_crlby > loc_uncer_max) | (loc_crlbz > loc_uncer_zmax) ;
zmask= loc_z > zmask_high | loc_z < zmask_low;
totmask = llmask | intmask | uncermask  | zmask;
totmask = zeros(size(llmask));
x_Real = (loc_x(~totmask))*x_pixelsize; % nm
y_Real = (loc_y(~totmask))*x_pixelsize; % nm
z_Real = (loc_z(~totmask))*z_pixelsize; % nm



% 6.3 VISPR Insitu PSF Reconstruction =========================================================================================================================================================================
for seq = 1:it_all
    if seq == it_all
        range = batchsize*(seq-1)+1:total_candidates;
    else
        range = batchsize*(seq-1)+1:batchsize*seq;
    end
    tic
    [t_Pcspline,CRLB,LogL]=mleFit_LM(subregion1(:, :, range), 6, 50,single(coeff_VISPR),sCMOSvarmap,1);  % 
    toc
%     Pcspline(range, :) =  [t_Pcspline,CRLB,LogL];
    Pcspline_XX(range, :) =  t_Pcspline;
end

result_VISPR = zeros(num_subregion1,5); % N, x, y, z, bg
result_VISPR(:,1) = Pcspline_XX(:,3);
result_VISPR(:,2) = Pcspline_XX(:,1);  % x
result_VISPR(:,3) = Pcspline_XX(:,2);  % y
result_VISPR(:,4) = Pcspline_XX(:,5);  % z
result_VISPR(:,5) = Pcspline_XX(:,4);
result_VISPR = single(result_VISPR);

x_VISPR = (x_start + result_VISPR(:,2) + delta_x_VISPR)*x_pixelsize; % pixel
y_VISPR = (y_start + result_VISPR(:,3) + delta_y_VISPR)*x_pixelsize; % pixel
z_VISPR = (result_VISPR(:,4)- size(coeff_VISPR,3)/2  + delta_z_VISPR)*z_pixelsize; % pixel

%% 7. show results
close all
x_GT = (x_final +(r_ROI+dis))*x_pixelsize;
y_GT = (y_final +(r_ROI+dis))*x_pixelsize;
% x_GT = (x_final)*x_pixelsize;
% y_GT = (y_final)*x_pixelsize;
z_GT = (z_final)*z_pixelsize;
% xshift_HF = mean(x_GT)-mean(x_HF);

show_pixel = 10;  % The pixel size of the final displayed super-resolution image

image_VISPR= show_full(show_pixel, z_VISPR,x_VISPR);
figure;
imshow(image_VISPR.^0.3,[]);title('VISPR PSF xz')

image_Real= show_full(show_pixel, z_Real,x_Real);
figure;
imshow(image_Real.^0.3,[]);title('Real PSF xz')

image_GTxz= show_full(show_pixel, z_GT,x_GT);
figure;
imshow(image_GTxz.^0.3,[]);title('Ground truth')

figure;
scatter(x_Real,z_Real,3,'filled');hold on
scatter(x_GT,z_GT,3,'filled');
legend('Real PSF','GT')
title('Real PSF vs GT')
axis equal
xlim([min(x_GT)-200,max(x_GT)+200])
ylim([min(z_GT)-200,max(z_GT)+200])

figure;
scatter(x_VISPR,z_VISPR,3,'filled');hold on
scatter(x_GT,z_GT,3,'filled');
legend('VISPR','GT')
title('VISPR vs GT')
axis equal
xlim([min(x_GT)-200,max(x_GT)+200])
ylim([min(z_GT)-200,max(z_GT)+200])

%% Store CSV results
frame = seg_display1.allcds_mask(:,3);
name = 'Reconstruction_result/EPFL';

name_GT = [name,'_GT'];
col = {'x[nm]','y [nm]','z [nm]'};
result_table = table(x_GT, y_GT, z_GT, 'VariableNames', col);
writetable(result_table, [name_GT,'.csv'])

name_Real_PSF = [name,'_Real_PSF'];
col = {'frame','x[nm]','y [nm]','z [nm]'};
result_table = table(frame, x_Real, y_Real, z_Real, 'VariableNames', col);
writetable(result_table, [name_Real_PSF,'.csv'])

VISPR_name_INSPR = [name,'_VISPR'];
col = {'frame','x[nm]','y [nm]','z [nm]'};
result_table = table(frame, x_VISPR, y_VISPR, z_VISPR, 'VariableNames', col);
writetable(result_table, [VISPR_name_INSPR,'.csv'])


id_GT = ones(numel(x_GT),1);
id_Real = ones(numel(x_Real),1).*2;
id_VISPR = ones(numel(x_VISPR),1).*3;

id_merge = cat(1,id_GT,id_Real,id_VISPR);
x_merge = cat(1,x_GT,x_Real,x_VISPR);
y_merge = cat(1,y_GT,y_Real,y_VISPR);
z_merge = cat(1,z_GT,z_Real,z_VISPR);
frame_GT = ones(numel(x_GT),1);
frame_merge = cat(1,frame_GT,frame,frame);
merge_name_INSPR = [name,'_merge']; 
col = {'id','frame','x[nm]','y [nm]','z [nm]'};
result_table = table(id_merge, frame_merge, x_merge, -1*y_merge, z_merge, 'VariableNames', col);
writetable(result_table, [merge_name_INSPR,'.csv'])