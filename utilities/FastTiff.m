function data = FastTiff(filename)
    warning('off','all') % Suppress all the tiff warnings
    tstack  = Tiff(filename);
    [I,J] = size(tstack.read());
    K = length(imfinfo(filename));
    data = zeros(I,J,K);
    data(:,:,1)  = tstack.read();
    for n = 2:K
        tstack.nextDirectory()
        data(:,:,n) = tstack.read();
    end
    warning('on','all')
end
