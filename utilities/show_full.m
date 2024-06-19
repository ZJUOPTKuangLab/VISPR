function img_full = show_full(pixelsize,xcoor, ycoor)

    % input的参数：xcoor和ycoor都是原始坐标
    tmpx = xcoor/pixelsize;
    tmpy = ycoor/pixelsize;

    tempx = round(tmpx - min(tmpx))+1;
    tempy = round(tmpy- min(tmpy))+1;
    img_full = zeros(max(tempx)+1, max(tempy)+1);
    for i = 1:size(tempx,1)
        img_full(tempx(i), tempy(i)) = img_full(tempx(i), tempy(i)) + 1; 
    end
    % figure;imshow(img_full,[]);title('img full')

end