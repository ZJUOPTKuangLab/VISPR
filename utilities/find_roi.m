function maxima = find_roi(PSF_data, thresh_dist, boxsz)

fs=2;        % filtersize
h = fspecial('gaussian',2*round(fs*3/2)+1,fs);


% mim = max(PSF_data,[],3);  % 找到中间位置
mim = PSF_data(:,:,round(size(PSF_data,3)/2));
mim = filter2(h,mim);
maxima = maximumfindcall(mim);

int = maxima(:,3);
mimc = mim(boxsz:end-boxsz,boxsz:end-boxsz);
bg = min(mimc(:));
if size(find(int > 2.5*bg), 1) >= 1
    cutoff = 2.5*bg;
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

end