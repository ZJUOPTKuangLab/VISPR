close all;
clear
clc
% We choose different zernike number to experiment. In most cases, we
% choose 21, 32, 45, 60 items for comparison. We found that the more items
% there are, the better the fitting effect is, but is more sensitive to noise. 
% Therefore, in general experiments, we choose 32 items as the default value.
global zernikenum;
zernikenum = 21;
addpath('source\')
addpath('utilities\')
addpath('shared\')

%% 1. Load SMLM raw data
ims = FastTiff('SMLM raw data\NPC-example.tif');
imageslicer(permute(ims(:,:,:),[2,1,3]));title('raw data')

zrange = -1200:10:1200; %z range of generated PSF
setup.is_imgsz = 1 ;
setup.is_sCMOS = 0;
setup.offset = 398.6;
setup.gain = 0.05;
setup.RefractiveIndex = 1.516;
setup.nMed = 1.335;
setup.Lambda = 680/1000;
setup.NA = 1.43;
setup.is_bg = 0;


%% 2. segmentation
boxsz = 19;
thresh_dist = 15;
thresh_bg = 5; % threshold needs to be carefully selected, otherwise the segmentation effect will not be ideal
setup.use_default_thresh = 1;
[subregion_ch1, seg_display] = crop_subregion_VISPR(ims, boxsz,thresh_dist,thresh_bg,setup);
num_subregion = size(subregion_ch1,3);
disp(['Finish segmentation! There are ' num2str(num_subregion) ' subregion images']);
subregion_ch1 = subregion_ch1(:,:,randi(num_subregion,1,min([num_subregion,5000])));

% image pre-processing
tmp = zeros(boxsz,boxsz);
for i = 1:size(subregion_ch1,3)
    tmp = subregion_ch1(:,:,i);
    BG = min(tmp(:));
    tmp = tmp - BG;
    subregion_ch1(:,:,i) = tmp;
end

imageslicer(single(subregion_ch1(:,:,:))); title('segmented PSF library')
mmed = quantile(subregion_ch1(:),0.3);
imbg = subregion_ch1(subregion_ch1<mmed);
bg =  mean(imbg(:));
Nph = mean(sum(sum(subregion_ch1-bg)));
Npixels = size(subregion_ch1,1);

%% 3. VISPR: generate in-situ PSF
zemit = -1000:50:1000;
iteration = 6;
empupil.min_similarity = 0.8;
empupil.bin_lowerBound = 50;
standard = 'Wyant';
run VISPR.m

parameters_insitu = set_parameters(31, numel(zrange), zrange, standard);
parameters_insitu.aberrations(:,3) = P(1:zernikenum);
insitu_PSF = gen_vectorial_psf(parameters_insitu);

imageslicer(insitu_PSF(:,:,1:10:end)); title('insitu PSF (interval:100nm)')
%% 4. Save spline coefficient of in-situ PSF for SMAP[1] localization
% 1. Ries, J. SMAP: a modular super-resolution microscopy analysis platform for SMLM data. Nat Methods 17, 870–872 (2020). https://doi.org/10.1038/s41592-020-0938-1
for i = 1:size(insitu_PSF,3) % 进行归一化
    tmp1 = insitu_PSF(:,:,i);
    tmp1 = tmp1/sum(tmp1(:));
    insitu_PSF(:,:,i) = tmp1;
end
coeff_VISPR = Spline3D_interp(insitu_PSF);
save('Cspline coeff\coeff_VISPR.mat','coeff_VISPR','insitu_PSF');




