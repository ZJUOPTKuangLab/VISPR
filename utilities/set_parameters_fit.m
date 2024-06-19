function paraFit = set_parameters_fit(PSF_stack, zemit, standard)
global zernikenum
paraFit.NA = 1.43;                                                    % numerical aperture of objective lens             
paraFit.refmed = 1.335;                                            % refractive index of sample medium
paraFit.refcov = 1.516;                                            % refractive index of converslip
paraFit.refimm = 1.516;                                            % refractive index of immersion oil
paraFit.NA = min(paraFit.NA, paraFit.refmed);
paraFit.lambda = 680;                                              % wavelength of emission
paraFit.pixelSizeX = 127;                                        % nm, pixel size of the image
paraFit.pixelSizeY = 116;                                        % nm, pixel size of the image

paraFit.zemit0 = 0;                                            % reference emitter z position, nm, distance of molecule to coverslip
paraFit.objStage0 = 0;                                            %  nm, initial objStage0 position,relative to focus at coverslip
paraFit.deltaz = zemit(1,end) - zemit(1,end-1);
paraFit.Npupil =64;                                             % sampling at the pupil plane

Npixels = size(PSF_stack,1); 
Nz = size(PSF_stack,3);
paraFit.sizeX = Npixels;                                        % number of pixels of the psf image
paraFit.sizeY = Npixels;                                        % number of pixels of the psf image
paraFit.sizeZ = Nz;

%initializing zernike polynomial coefficients
% standard = 'Wyant';  % 'Wyant', 'OSA'
paraFit.aberrations = Select_Zernike_terms(zernikenum, standard);
numAberrations = size(paraFit.aberrations,1);
shared = [ones(1,numAberrations) 1 1 1 0 0];  % 1 is shared parameters between z slices, 0 is free parameters between z slices, only consider  [x, y, z, I, bg]

sumShared = sum(shared);
numparams = (zernikenum+5) * Nz-sumShared*(Nz-1);
thetainit = zeros(numparams,1);
bg = zeros(1,Nz);
Nph = zeros(1,Nz);
x0 = zeros(1,Nz);
y0 = zeros(1,Nz);
z0 = zeros(1,Nz);

% center of mass with nm unit
ImageSizex = paraFit.pixelSizeX*Npixels/2;
ImageSizey = paraFit.pixelSizeY*Npixels/2;

DxImage = 2*ImageSizex/paraFit.sizeX;
DyImage = 2*ImageSizey/paraFit.sizeY;
ximagelin = -ImageSizex+DxImage/2:DxImage:ImageSizex;
yimagelin = -ImageSizey+DyImage/2:DyImage:ImageSizey;
[YImage,XImage] = meshgrid(yimagelin,ximagelin);

for i = 1:Nz
    dTemp = PSF_stack(:,:,i);
    mmed = quantile(dTemp(:),0.7);
    imbg = dTemp(dTemp<mmed);
    bg(i) = mean(imbg(:));
    bg(i) = max(bg(i),1);
    Nphim = dTemp - bg(i);
    Nphim(Nphim<0) = 1e-6;

    Nph(i) = sum(sum(Nphim));
    x0(i) = sum(sum(XImage.*dTemp))/Nph(i);
    y0(i) = sum(sum(YImage.*dTemp))/Nph(i);
    z0(i) = 0;%
end

allTheta = zeros(numAberrations+5,Nz);
allTheta(numAberrations+1,:)=x0';
allTheta(numAberrations+2,:)=y0';
allTheta(numAberrations+3,:)=z0';
allTheta(numAberrations+4,:)=Nph';
allTheta(numAberrations+5,:)=bg';
allTheta(1:numAberrations,:) = repmat(paraFit.aberrations(:,3),[1 Nz]);

map = zeros(numparams,3);
n=1;
for i = 1:numAberrations+5
    if shared(i)==1
        map(n,1)= 1;
        map(n,2)=i;
        map(n,3)=0;
        n = n+1;
    elseif shared(i)==0
        for j = 1:Nz
            map(n,1)=0;
            map(n,2)=i;
            map(n,3)=j;
            n = n+1;
        end
    end
end


for i = 1:numparams
    if map(i,1)==1
        thetainit(i)= mean(allTheta(map(i,2),:));
    elseif map(i,1)==0
         thetainit(i) = allTheta(map(i,2),map(i,3));
    end
end

% we assume that parameters for zernike coefficiens are always linked
zernikecoefsmax = 0.25*paraFit.lambda*ones(numAberrations,1);
paraFit.maxJump = [zernikecoefsmax',paraFit.pixelSizeX*ones(1,max(Nz*double(shared(numAberrations+1)==0),1)),paraFit.pixelSizeY*ones(1,max(Nz*double(shared(numAberrations+2)==0),1)),500*ones(1,max(Nz*double(shared(numAberrations+3)==0),1)),2*max(Nph(:)).*ones(1,max(Nz*double(shared(numAberrations+4)==0),1)),100*ones(1,max(Nz*double(shared(numAberrations+5)==0),1))];
paraFit.numparams = numparams;
paraFit.numAberrations = numAberrations;

paraFit.zemitStack = zemit; % move emitter
tempGauss1 = 0;

paraFit.objStageStack = zeros(size(PSF_stack,3),1);
paraFit.ztype = 'emitter';
paraFit.map = map;
paraFit.Nitermax = 75;

paraFit.PSF_stack = PSF_stack;
paraFit.thetainit = thetainit;
paraFit.shared = shared;
paraFit.tempGauss1 = tempGauss1;

% Bead parameters for convolution with PSF and derivatives, beaddiameter in nm
paraFit.bead = true;
% parameters.beaddiameter = 175;
paraFit.beaddiameter = 100;

% check on meaningfullness of bead convolution
if paraFit.beaddiameter<=paraFit.pixelSizeX
  paraFit.bead = false;
end


end