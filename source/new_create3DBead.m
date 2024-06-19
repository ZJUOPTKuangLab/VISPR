function bead = new_create3DBead(parameters,bead_xsize,bead_zsize,zoom)


x_pixelsize = parameters.pixelSizeX;
y_pixelsize = parameters.pixelSizeY;
z_pixelsize = parameters.deltaz;
D = parameters.beaddiameter;

bead = get_realbead(bead_xsize, bead_zsize, x_pixelsize, y_pixelsize, z_pixelsize, zoom, D);
real_xsize = parameters.sizeX;
real_zsize = parameters.sizeZ;
padding_xsize = (real_xsize-bead_xsize)/2;
padding_zsize = floor((real_zsize-bead_zsize)/2);
bead = padarray(bead,[padding_xsize padding_xsize padding_zsize],0,'both');

end