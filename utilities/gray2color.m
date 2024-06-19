function color_map = gray2color(heat_map)
%将灰度图像转换为彩色图像
cmap = turbo;
sz  = size(heat_map);
color_map = zeros(sz(1), sz(2), 3);
for i = 1:256
    ind = find(heat_map == i-1);
    [r,c] = ind2sub(sz, ind);
    for j = 1:3
        color_map(sub2ind([sz(1), sz(2), 3], r, c, j+zeros(size(r)))) = cmap(i, j);
    end 
end
end

