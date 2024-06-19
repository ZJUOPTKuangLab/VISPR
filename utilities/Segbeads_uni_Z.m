function [PSF_stack, maxima] = Segbeads_uni_Z(PSF_data, boxsz, thresh_dist)

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
if size(find(int > 5*bg), 1) >= 1
    cutoff = 5*bg;
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

index_outrange = [];
for l=1:size(maxima,1)
    pos=maxima(l,1:2);
    if pos(2)+roisizeh < size(imstack,1) && pos(1)+ roisizeh < size(imstack,2)...
            && pos(2)-roisizeh > 0 && pos(1)-roisizeh > 0
    else
       index_outrange = [index_outrange,l];
    end
end
maxima(index_outrange,:) = [];


PSF_stack_dataset = cell(size(maxima,1),1);
for l=1:size(maxima,1)
    pos=maxima(l,1:2);
    if pos(2)+roisizeh < size(imstack,1) && pos(1)+ roisizeh < size(imstack,2)...
            && pos(2)-roisizeh > 0 && pos(1)-roisizeh > 0
       PSF_stack = PSF_stack + imstack(pos(2)+rsr,pos(1)+rsr,:);
       PSF_stack_dataset{l,1} = imstack(pos(2)+rsr,pos(1)+rsr,:);    
    end
end



for l=1:size(maxima,1)
    PSF_image = PSF_stack_dataset{l,1};
    if isempty(PSF_image) == 0
        imageslicer(permute(PSF_image,[2,1,3]))
        title(['spot: ', num2str(l)]);
    end
end

PSF_stack = PSF_stack./size(maxima,1);


end