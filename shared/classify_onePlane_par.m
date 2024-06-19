%% This function classifies single molecule to giving reference 
%  Z-position images, and average these Z-postion images in a certain
%  group.
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, October 2019
% 
function [ims_Zcal_ave_plane1, index_record_Zplanes] = classify_onePlane_par(subregion_ch1, ref_plane1, empupil) 
% 注：此处的ref_plane1是其实是乘了5000的接近归一化的图像(因为还乘上了OTF rescale)

disp('Calculate the similarity between reference images and single molecules');
img_plane1 = single(subregion_ch1);    %normalization image

num_img = size(img_plane1,3);  %3535个候选点
num_ref = size(ref_plane1, 3);   % 21个reference z position
imsz = empupil.imsz;    % image size of the ROI psf 

similarity_in_Plane1 = zeros(num_ref, num_img);  %similarity  21*3535

index_similarity = zeros(num_ref, num_img); %index
shift_row_plane1 = zeros(num_ref, num_img); %X,Y shift in two planes
shift_col_plane1 = zeros(num_ref, num_img);


parfor ii = 1 : num_img     % 2023年5月13日，这里最开始是parfor 为了调试代码

    for jj = 1 : num_ref
        [shift1, shift2, tmpval1] = registration_in_each_channel(ref_plane1(:,:,jj), img_plane1(:,:,ii));
        % 注意：在找到similarity tmpval1的时候，有对ref image进行归一化(减去均值，并且除以std，具体见该函数的最后一行cc2.m function)，所以其实前面生成的scalePSF没有能量归一化问题不大
        % 并且，在计算img同ref_img的similarity的时候，其实是先生成img_shift(即是以ref_img的中心为中心)进行XY
        % registration，然后再计算相关的similarity

        if (abs(shift1) > 6 || abs(shift2) > 6)  %this value can be changed, now fixed
            continue;  % continue的意思是直接进行下一步循环了
        end     % 所以这一步是filter掉大于一定pixelshift的数据，但是可以更改
        
        %get X, Y shift in plane1
        shift_row_plane1(jj, ii) = shift1;
        shift_col_plane1(jj, ii) = shift2;
            
        if tmpval1 < empupil.min_similarity
            continue;
        end  % 这一步是filter out掉小于最小similarity的数据
     
        similarity_in_Plane1(jj, ii) = tmpval1;
    end    
end


%%
for ii = 1 : num_img
    % Sort the similarity
    [sort_similarity, index_sort] =  sort(similarity_in_Plane1(:,ii),'descend'); 
    if (sort_similarity(1) == 0)
        continue;    % 这一段的作用是，如果这个PSF每一个reference都不是的话，那么就直接进行到下一个点了
    end
    % Determine single molecule image belong which reference image
    index_similarity(index_sort(1), ii) = 1;   % index similarity也是21*3535的矩阵
    for jj = 2 : num_ref  % num_ref = 21
        if (sort_similarity(jj) >= sort_similarity(1)-0.0) && (abs(index_sort(jj)-index_sort(1)) == 1)    %now fixed
            index_similarity(index_sort(jj),ii) = 1;  %相当于是为了防止不同平面的similarity相同，并且是相邻两个平面的候选点
        else
            break;
        end
    end
end


%% Updata average images
disp('Update average images');

ims_Zcal_ave_plane1 = zeros(imsz,imsz,num_ref);   % 32*32*21的三维矩阵
index_record_Zplanes = zeros(num_ref,1);    % 21*1的矩阵，记录哪些位置的是有相应的平均值图像的
for ii = 1 : num_ref    % 当i=6的时候，才第一次开始sz_index>bin_lowerBound了
    index_selection = find(index_similarity(ii,:) == 1);
    sz_index = size(index_selection,2);
    if sz_index > empupil.bin_lowerBound    %Edited by FX, each group must have enough spots 
        ims_plane1_shift = zeros(imsz,imsz,sz_index);  %32*32*sz_index 
        for jj = 1 : sz_index
            ims_plane1_shift(:,:,jj) = FourierShift2D(similarity_in_Plane1(ii, index_selection(jj)) .* img_plane1(:,:,index_selection(jj)), [shift_row_plane1(ii, index_selection(jj)) shift_col_plane1(ii, index_selection(jj))]);            
        end   
        % 这一步是通过相似度的大小进行图像的加权，然后进行平均，但是我记得前面已经生成过了相应的平移图了
        % 并且这里的img_plane1就是 single(subregion_ch1)   %normalization image

        % average the images
        index_record_Zplanes(ii) = 1;
        ims_Zcal_ave_plane1(:,:,ii) = mean(ims_plane1_shift,3);       % 所以这里的average的均值基本应该就是0
    end
end


%% Remove too far away Reassembled PSF
%Note：第一个for循环是从中间位置(焦面)往负方向走，如果和距离焦面更近位置的间隔差的太大，那么就变为去掉
tmp_record_keep = find(index_record_Zplanes == 1);
size_tmp = length(tmp_record_keep);

for ii = ceil(size_tmp/2) : -1 : 2
    if tmp_record_keep(ii) >= tmp_record_keep(ii-1)+2   % 那就是至少差3
        index_record_Zplanes(tmp_record_keep(ii-1)) = 0;
        tmp_record_keep(ii-1) = tmp_record_keep(ii);
    end
end

for ii = ceil(size_tmp/2) : size_tmp-1
    if tmp_record_keep(ii) <= tmp_record_keep(ii+1)-2
        index_record_Zplanes(tmp_record_keep(ii+1)) = 0;
        tmp_record_keep(ii+1) = tmp_record_keep(ii);
    end
end 
