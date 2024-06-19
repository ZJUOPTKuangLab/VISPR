function [subregion_ch1,seg_display] = crop_subregion_VISPR(data1, boxsz, thresh_dist, thresh, setup)

if setup.is_sCMOS   %sCMOS case 
    % sCMOS parameters
    offsetim_ch1 = repmat(setup.sCMOS_input.ccdoffset_ch1,[1 1 size(data1,3)]);    
    gainim_ch1 = repmat(setup.sCMOS_input.gain_ch1,[1 1 size(data1,3)]);   
    data1_in = (data1 - offsetim_ch1) ./ gainim_ch1;
else    %EMCCD case
    data1_in = (data1 - setup.offset) /setup.gain;
%     data1_in = (data1 - setup.offset); 
end

data1_in(data1_in<=0) = 1e-6;


imsz_original = size(data1_in);
if setup.is_imgsz == 1
    rangemin=[5, 5];
    rangemax=[imsz_original(1)-5, imsz_original(2)-5];
else 
    rangemin=[1, 1];
    rangemax=[imsz_original(1), imsz_original(2)];
end

ims_ch1 = single(data1_in(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
% ims_ch1 = data1_in;
ims_detect = ims_ch1;


%% Detect single molecules
disp('obtain sub_region centers');
roisizeh = round((boxsz-1)/2); %create extra space if we need to shift;
rsr = -roisizeh:roisizeh;
fs=2;        % filtersize
h = fspecial('gaussian',2*round(fs*3/2)+1,fs);
allcds = [];
use_default_thresh = setup.use_default_thresh;
parfor i = 1:size(ims_detect,3)
    ims = ims_detect(:,:,i);
    mim = filter2(h,ims);
    
    % find the local maximum point
    maxima = maximumfindcall(mim);
    int = maxima(:,3);
    mimc = mim(boxsz:end-boxsz,boxsz:end-boxsz);
    bg = min(mimc(:));
    
    if size(find(int > thresh*bg), 1) >= 1 && use_default_thresh == 0
        cutoff = thresh*bg;
    else
        disp('use default thresh')
        mmed = quantile(mimc(:),0.9);
        imt = mimc(mimc<mmed);
        sm = sort(int);
        mv = mean(sm(end-5:end));
        cutoff = mean(imt(:)) + max(2.5*std(imt(:)),(mv-mean(imt(:)))/15);
    end
    
    if any(int>cutoff)
        maxima=maxima(int>cutoff,:);
    else
        [~,indm]=max(int);
        maxima=maxima(indm,:);
    end
    
    % Remove too close spots
    pos = maxima(:,1:2);
    psf_dist = pdist(pos);
    num_tmp = 1;
    psf_remove = [];
    for ii = 1 : (size(pos, 1)-1)
        for jj = ii + 1 : size(pos, 1)
            if (psf_dist(num_tmp) < thresh_dist)
                psf_remove = [psf_remove ii jj];
            end
            num_tmp = num_tmp + 1;
        end
    end
    maxima(psf_remove,:) = [];
    maxima(:,3) = [];

    A = ones(size(maxima,1),2);
    A(:,1) = i-1;

    maxima = [maxima, A];
    allcds = [allcds; maxima];
end
%% selected spots

imsize_x = size(ims_detect, 1);
imsize_y = size(ims_detect, 2);

boundmask=(allcds(:,1)<=boxsz/2)|(allcds(:,1)>=imsize_x-boxsz/2)|(allcds(:,2)<=boxsz/2)|(allcds(:,2)>=imsize_y-boxsz/2) | allcds(:,4) == 0;
allcds_mask = allcds(~boundmask,:);

%% Segmentation sub-region of single molecule

ims_ch1 = single(data1_in(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));

[subregion_ch1, l1, t1] = cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch1);



%% save display
seg_display.ims_ch1 = ims_ch1;
seg_display.allcds_mask = allcds_mask;
seg_display.t1 = t1;
seg_display.l1 = l1;
