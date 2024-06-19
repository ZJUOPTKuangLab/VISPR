function bead = get_realbead(xsize, zsize, x_pixelsize, y_pixelsize, z_pixelsize, zoom, D)
    bead = zeros(xsize, xsize, zsize);

    xsize_eq = xsize*zoom;
%     ysize_eq = xsize*zoom;
    zsize_eq = zsize*zoom;
    xmid = floor((xsize_eq+1)/2);
%     ymid= xmid;
    zmid = floor((zsize_eq+1)/2);
    xx = 1-xmid:xsize_eq-xmid;
    yy = xx;
    zz = 1-zmid:zsize_eq-zmid;

    [ys, xs, zs] = meshgrid(yy,xx,zz);
    xs = xs*x_pixelsize/zoom;
    ys = ys*y_pixelsize/zoom;
    zs = zs*z_pixelsize/zoom;

    R_beads = D/2;

    bead_eq = xs.^2+ys.^2+zs.^2<R_beads.^2;
    for i = 1:xsize
        for j = 1:xsize
            for k = 1:zsize
                tmp = bead_eq( (i-1)*zoom+1:i*zoom, (j-1)*zoom+1:j*zoom, (k-1)*zoom+1:k*zoom);
                bead(i, j, k) =sum(tmp(:));
            end
        end
    end

    bead = bead./sum(bead(:));

end