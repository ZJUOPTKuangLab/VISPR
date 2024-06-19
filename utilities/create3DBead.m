function bead = create3DBead(parameters)
% This function creates an image of a bead for a 3D convolution with a PSF.


zrange = parameters.zemit;
pixelSizeX = parameters.pixelSizeX;
pixelSizeY = parameters.pixelSizeY;
sizeX = parameters.sizeX;
sizeY = parameters.sizeY;
ImageSizex = pixelSizeX*sizeX/2; 
ImageSizey = pixelSizeY*sizeY/2; 

zmin = min(zrange);
zmax = max(zrange);
Mx = parameters.sizeX;
My = parameters.sizeY;
Mz = parameters.sizeZ;
beaddiameter = parameters.beaddiameter;

% image coordinate sampling (in physical length units).
DxImage = 2*ImageSizex/Mx;
DyImage = 2*ImageSizey/My;
ximagelin = -ImageSizex+DxImage/2:DxImage:ImageSizex;
yimagelin = -ImageSizey+DyImage/2:DyImage:ImageSizey;
if Mz==1
    zimagelin = 0;
else
zimagelin = linspace(zmin,zmax,Mz);
end
[YImage,XImage,ZImage] = meshgrid(yimagelin,ximagelin,zimagelin);

% compute coordinates in spherical coordinates 
rho = sqrt(XImage.^2+YImage.^2+ZImage.^2);

bead = rho<(beaddiameter/2);

bead = bead/sum(sum(sum(bead)));

end

