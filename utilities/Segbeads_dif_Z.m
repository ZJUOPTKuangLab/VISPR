function [PSF_stack, maxima] = Segbeads_dif_Z(PSF_data, boxsz, thresh_dist)

% Segmentation
roisizeh = round((boxsz-1)/2); % create extra space if we need to shift;
rsr = -roisizeh:roisizeh;
znum = size(PSF_data,3);
Npixels = boxsz;
PSF_stack = zeros(Npixels,Npixels,znum);

fs=2;        % filtersize
h = fspecial('gaussian',2*round(fs*3/2)+1,fs);

% Select the bright points in the image
imstack = PSF_data;
mim = max(imstack,[],3);
mim = filter2(h,mim);
maxima = maximumfindcall(mim);
int = maxima(:,3);
mimc = mim(boxsz:end-boxsz,boxsz:end-boxsz);
bg = min(mimc(:));
if size(find(int > 1.7*bg), 1) >= 1
    cutoff = 1.7*bg;
else
    mmed = quantile(mimc(:),0.8);
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


psf_num = size(maxima,1);
znum = size(PSF_data,3);

for i = 1:znum
    raw_bead = PSF_data(:,:,i);
    % bead = raw_bead-min(raw_bead(:)); %fast fix for offset;
    bead = raw_bead;

    mim = filter2(h,bead);
    maxima = maximumfindcall(mim);
    int = maxima(:,3);
    mimc = mim(boxsz:end-boxsz,boxsz:end-boxsz);
    bg = min(mimc(:));

    if size(find(int > 1.5*bg)) >=1
        cutoff = 1.5*bg;
    else
        mmed = quantile(mimc(:),0.8);
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

    maxima = sortrows(maxima,3,'descend');
    singlepsf = zeros(boxsz,boxsz);

    maxima = maxima(1:psf_num,:);

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

    for l=1:size(maxima,1)
        pos=maxima(l,1:2);
        if pos(2)+roisizeh < size(bead,1) && pos(1)+ roisizeh < size(bead,2)...
            && pos(2)-roisizeh > 0 && pos(1)-roisizeh > 0

            singlepsf = singlepsf + bead(pos(2)+rsr,pos(1)+rsr);
        end
    end
    psf = singlepsf./size(maxima,1);
    PSF_stack(:,:,i) = psf;

end



end