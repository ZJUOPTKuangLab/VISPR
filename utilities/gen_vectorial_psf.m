function [PSF] = gen_vectorial_psf(parameters)
% This function calculates the pupil matrix Q_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component. Then it sums the six components of pupil matrix to get the
% vectoorial psf model based on parameters.

% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% wavelength (in nm), sampling in pupil
NA = parameters.NA;
refmed = parameters.refmed;
refcov = parameters.refcov;
refimm = parameters.refimm;
lambda = parameters.lambda;
Npupil = parameters.Npupil;  % sampling at the pupil plane
% 计算图像pixel的一些参数
sizeX = parameters.sizeX;    % number of pixels of the PSF image in X
sizeY = parameters.sizeY;    % number of pixels of the PSF image in Y
sizeZ = parameters.sizeZ;    % Z轴stack数，比如21

% center of mass with nm unit
ImageSizex = (parameters.pixelSizeX*parameters.sizeX/2);
ImageSizey = (parameters.pixelSizeY*parameters.sizeY/2);
DxImage = 2*ImageSizex/sizeX;
DyImage = 2*ImageSizey/sizeY;


% enlarging ROI in order to accommodate convolution with bead in
% computation PSF and PSFderivatives
if isfield(parameters,'bead')
  if parameters.bead == true
    beaddiameter = parameters.beaddiameter;
    DeltaMx = 2*ceil((beaddiameter/2)/DxImage);
    sizeX = sizeX + DeltaMx;
    ImageSizex = ImageSizex+DeltaMx*DxImage/2;
    DeltaMy = 2*ceil((beaddiameter/2)/DyImage);
    sizeY = sizeY + DeltaMy;
    ImageSizey = ImageSizey+DeltaMy*DyImage/2;
  end
end
ImageSizex = ImageSizex*NA/lambda;
ImageSizey = ImageSizey*NA/lambda;
zemit = parameters.zemit;
zmin = min(zemit);
zmax = max(zemit);
DzImage = (zmax - zmin)/(sizeZ-1); 

% enlarging axial range in order to accommodate convolution with bead
if isfield(parameters,'bead')
    if parameters.bead == true
        beaddiameter = parameters.beaddiameter;
        DeltaMz = 2*ceil((beaddiameter/2)/DzImage);
        sizeZ = sizeZ + DeltaMz;
        zmin = min(zemit)-DeltaMz*DzImage/2;
        zmax = max(zemit)+DeltaMz*DzImage/2;
    end
    zemit = zmin:DzImage:zmax;
end

% pupil radius (in diffraction units) and pupil coordinate sampling
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,sizeX);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,sizeY);

% calculation of relevant Fresnel-coefficients for the interfaces
% between the medium and the cover slip and between the cover slip
% and the immersion fluid
CosThetaMed  = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
CosThetaCov = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2);
CosThetaImm = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2);
%need check again, compared to original equation, FresnelPmedcov is
%multipiled by refmed
FresnelPmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaCov+refcov*CosThetaMed);
FresnelSmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaMed+refcov*CosThetaCov);
FresnelPcovimm = 2*refcov*CosThetaCov./(refcov*CosThetaImm+refimm*CosThetaCov);
FresnelScovimm = 2*refcov*CosThetaCov./(refcov*CosThetaCov+refimm*CosThetaImm);
FresnelP = FresnelPmedcov.*FresnelPcovimm;
FresnelS = FresnelSmedcov.*FresnelScovimm;

% Apoidization for sine condition
% apoid = sqrt(CosThetaImm)./CosThetaMed;
apoid = 1 ./ sqrt(CosThetaMed);
% definition aperture
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
Amplitude = ApertureMask.*apoid;

% setting of vectorial functions
Phi = atan2(YPupil,XPupil);  %zernike多项式的Φ变量,也就是方位角
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
SinTheta = sqrt(1-CosTheta.^2);

%P,S偏振方向对于x，y，z方向的六个分量
pvec{1} = FresnelP.*CosTheta.*CosPhi;
pvec{2} = FresnelP.*CosTheta.*SinPhi;
pvec{3} = -FresnelP.*SinTheta;
svec{1} = -FresnelS.*SinPhi;
svec{2} = FresnelS.*CosPhi;
svec{3} = 0;
PolarizationVector = cell(2,3);
for jtel = 1:3
  PolarizationVector{1,jtel} = CosPhi.*pvec{jtel}-SinPhi.*svec{jtel};
  PolarizationVector{2,jtel} = SinPhi.*pvec{jtel}+CosPhi.*svec{jtel};
end

% calculate wavevector inside immersion fluid and z-component inside medium
% 计算波矢
wavevector = cell(1,3);
wavevector{1} = (2*pi*NA/lambda)*XPupil;
wavevector{2} = (2*pi*NA/lambda)*YPupil;
wavevector{3} = (2*pi*refimm/lambda)*CosThetaImm;
wavevectorzmed = (2*pi*refmed/lambda)*CosThetaMed; 

% calculation aberration function native to optical system
% 使用Zernike多项式表达，Zernike系数单位是nm
Waberration = zeros(size(XPupil));
orders = parameters.aberrations(:,1:2);
zernikecoefs = squeeze(parameters.aberrations(:,3));
normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0))); %归一化参数
zernikecoefs = normfac.*zernikecoefs;
allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
for j = 1:numel(zernikecoefs)
  Waberration = Waberration+zernikecoefs(j)*squeeze(allzernikes(j,:,:));  
end
% compute effect of refractive index mismatch
% 除了系统的Zernike带来的相位偏差外还要包含样品，油，盖片的折射率不匹配带来的相位偏差
% 原理还没彻底弄懂，就先不写
% 相位因子算出来
PhaseFactor = exp(2*pi*1i*Waberration/lambda);  %除以lamdazai再乘上2pi是为了把单位从波长转成rad


% calculate intensity normalization function using the PSFs at focus
% position without any aberration 用于归一化
FieldMatrix = cell(2,3);
for itel = 1:2
  for jtel = 1:3
    PupilFunction = Amplitude.*PolarizationVector{itel,jtel};
    IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
    FieldMatrix{itel,jtel} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
  end
end

intFocus = zeros(sizeX,sizeY);
for jtel = 1:3
    for itel = 1:2
        intFocus = intFocus + (1/3)*abs(FieldMatrix{itel,jtel}).^2;
    end
end
normIntensity = sum(intFocus(:));

%% FieldMatrix/FieldMatrixDerivatives
xemit = ones(1,sizeZ)* parameters.xemit;
yemit = ones(1,sizeZ)* parameters.yemit;

objStage0= parameters.objStage0;  % nm, distance from coverslip to nominal focal plane, when RI match,this is useless
zemit0 = parameters.zemit0;  % reference emitter z position, nm, distance of molecule to coverslip

zemitStack = zeros(sizeZ,1)';  % move emitter
objStageStack = zeros(sizeZ,1)'; % move objStage
ztype = parameters.ztype;

numAberrations = size(parameters.aberrations,1);
numders = 3+numAberrations;

FieldMatrix = cell(2,3,sizeZ);
FieldMatrixDerivatives = cell(2,3,sizeZ,numders);

% dudxyzdzer_field = cell(2,3,sizeZ,3,numAberrations);


if strcmp(ztype,'stage')
    objStage = zemit + objStageStack(:)+objStage0(:);
    zemit = zemitStack(:)+zemit0(:);
elseif strcmp(ztype,'emitter')
    objStage = objStageStack(:) + objStage0(:);
    zemit = zemit' + zemitStack(:)+ zemit0(:);
end


for jz = 1:sizeZ
    % phase contribution due to position(xyz) of the emitter；计算xy位置导致的相位
    Wxyz= (-1*xemit(jz))*wavevector{1}+(-1*yemit(jz))*wavevector{2}+(zemit(jz))*wavevectorzmed;
    PositionPhaseMask = exp(1i*(Wxyz+(objStage(jz))*wavevector{3}));

    for itel = 1:2
        for jtel = 1:3
            % pupil functions and FT to matrix elements
            PupilMatrix = Amplitude.*PhaseFactor.*PolarizationVector{itel,jtel};
            PupilFunction = PositionPhaseMask.*PupilMatrix;
            IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
            FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

            % pupil functions for xy-derivatives and FT to matrix elements
            for jder = 1:2
                PupilFunction_xy = -1i*wavevector{jder}.*PositionPhaseMask.*PupilMatrix;
                IntermediateImage = transpose(cztfunc(PupilFunction_xy,Ay,By,Dy));
                FieldMatrixDerivatives{itel,jtel,jz,numAberrations+jder} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
            end

            % pupil functions for z-derivative and FT to matrix elements
            if strcmp(ztype,'stage') % find stage position，用的imm的矢量
                PupilFunction_z = 1i*wavevector{3}.*PositionPhaseMask.*PupilMatrix;
            elseif strcmp(ztype,'emitter') % find emitter position，用的是med的矢量
                PupilFunction_z = 1i*wavevectorzmed.*PositionPhaseMask.*PupilMatrix;
            end
            IntermediateImage = transpose(cztfunc(PupilFunction_z,Ay,By,Dy));
            FieldMatrixDerivatives{itel,jtel,jz,numAberrations+3} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

            % pupil functions for Zernike mode-derivative and FT to matrix elements
            %        n = 1;
            for jzer = 1:numAberrations
                PupilFunction_zer = (2*pi*1i*normfac(jzer)*squeeze(allzernikes(jzer,:,:))/lambda).*PositionPhaseMask.*PupilMatrix;
                IntermediateImage = transpose(cztfunc(PupilFunction_zer,Ay,By,Dy));
                FieldMatrixDerivatives{itel,jtel,jz,jzer} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
            end
           
            % dfield/(dxyzdzer)
%             for jpos = 1:3
%                 if jpos~=3
%                     for jzer = 1:numAberrations
%                         tmp = wavevector{jpos}.*squeeze(allzernikes(jzer,:,:))*2*pi*normfac(jzer)/lambda.*PositionPhaseMask.*PupilMatrix;
%                         IntermediateImage = transpose(cztfunc(tmp,Ay,By,Dy));
%                         dudxyzdzer_field{itel,jtel,jz,jpos,jzer}=transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
%                     end
%                 else
%                     for jzer = 1:numAberrations
%                         if strcmp(ztype,'stage')
%                             tmp = -1*wavevector{3}.*squeeze(allzernikes(jzer,:,:))*2*pi*normfac(jzer)/lambda.*PositionPhaseMask.*PupilMatrix;
%                         elseif strcmp(ztype,'emitter') 
%                             tmp = -1*wavevectorzmed.*squeeze(allzernikes(jzer,:,:))*2*pi*normfac(jzer)/lambda.*PositionPhaseMask.*PupilMatrix;
%                         end
%                         IntermediateImage = transpose(cztfunc(tmp,Ay,By,Dy));
%                         dudxyzdzer_field{itel,jtel,jz,jpos,jzer}=transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
%                     end                    
%                 end
%             end

        end
    end
end

%% get 3D vectorial PSF
PSF = zeros(sizeX,sizeY,sizeZ);
PSFderivatives = zeros(sizeX,sizeY,sizeZ,numders);
% dudxyzdzer = zeros(sizeX,sizeY,sizeZ,3,numAberrations);

for jz = 1:sizeZ
    for jtel = 1:3
        for itel = 1:2
            PSF(:,:,jz) = PSF(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
            for jder = 1:numders
                PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) +...
                    (2/3)*real(conj(FieldMatrix{itel,jtel,jz}).*...
                    FieldMatrixDerivatives{itel,jtel,jz,jder});
            end
            %% du/(dxyzdzer)
%             for jpos=1:3
%                 for jzer = 1:numAberrations
%                     dudxyzdzer(:,:,jz,jpos,jzer) = dudxyzdzer(:,:,jz,jpos,jzer)+...
%                         (2/3)*real( conj(FieldMatrixDerivatives{itel,jtel,jz,jzer}).*...
%                         FieldMatrixDerivatives{itel,jtel,jz,numAberrations+jpos}+...
%                         conj(FieldMatrix{itel,jtel,jz}).*...
%                         dudxyzdzer_field{itel,jtel,jz,jpos,jzer});                  
%                 end
%             end
            %%
        end
    end
end

PSF = PSF/normIntensity;
PSFderivatives = PSFderivatives/normIntensity;
% dudxyzdzer =dudxyzdzer / normIntensity;

% 3D convolution of the PSFs and derivatives with a bead
if isfield(parameters,'bead')
  if parameters.bead == true
    bead = new_create3DBead(parameters,5,5,101);
    PSF = convn(bead,PSF,'same');
    [Mx,My,Mz] = size(PSF);
    tempderivs = zeros(Mx,My,Mz,numders);
    for jder = 1:numders
      tempderivs(:,:,:,jder) = convn(bead,squeeze(PSFderivatives(:,:,:,jder)),'same');
    end
    PSFderivatives = tempderivs;
  end
end


end











