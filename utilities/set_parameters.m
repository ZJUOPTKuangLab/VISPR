function parameters = set_parameters(Npixels, Nz, zemit, standard)
global zernikenum
parameters.NA = 1.43;
parameters.refmed = 1.335;
parameters.refcov = 1.516;
parameters.refimm = 1.516;
parameters.NA = min(parameters.NA, parameters.refmed);
parameters.lambda = 680;
parameters.pixelSizeX = 127;       % nm, pixel size of the image
parameters.pixelSizeY = 116;       % nm, pixel size of the image

parameters.sizeZ = Nz;              % number of discrete z positions
parameters.objStage0 = 0;              % nm, distance from coverslip to nominal focal plane, when RI match,this is useless
parameters.zemit0 = -1*parameters.refmed/parameters.refimm*(parameters.objStage0);% reference emitter z position, nm, distance of molecule to coverslip
parameters.zemitStack = zeros(Nz,1)';  % move emitter,在这个函数中不使用
parameters.objStageStack = zeros(Nz,1)'; %move objStage,在这个函数中不使用
parameters.ztype = 'emitter';
parameters.deltaz = zemit(1,end) - zemit(1,end-1);
                        
parameters.sizeX = Npixels;             % number of pixels of the psf image
parameters.sizeY = Npixels;             % number of pixels of the psf image
parameters.Npupil = 64;                 % sampling at the pupil plane

% total 21 zernike aberration orders, n,m,amplitude(nm)
% standard = 'Wyant'; % Wyant, OSA
parameters.aberrations = Select_Zernike_terms(zernikenum, standard);
parameters.aberrations(1,3) = 80;  % Vertical Astig
parameters.aberrations(2,3) = 0; % Oblique Astig


parameters.xemit = 0;                %nm x positions of each PSF
parameters.yemit = 0;                %nm y positions of each PSF
parameters.zemit = zemit;                         %nm z positions of each PSF
parameters.objStage = linspace(-1000,1000,Nz)*0;   %nm objective positions of each PSF

Nphotons = 30000 +5000*rand(1,Nz);    % photons for each PSF
bg = 100 + 100*rand(1,Nz);             % background for each PSF
parameters.Nphotons = Nphotons;
parameters.bg = bg;

% Bead parameters for convolution with PSF and derivatives, beaddiameter in nm
parameters.bead = true;
% parameters.beaddiameter = 175;
parameters.beaddiameter = 100;

% check on meaningfullness of bead convolution
if parameters.beaddiameter<=parameters.pixelSizeX
  parameters.bead = false;
end

end