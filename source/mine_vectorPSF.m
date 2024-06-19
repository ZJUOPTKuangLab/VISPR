function [dudt,model] =  mine_vectorPSF(theta,parameters,map,shared,bead)
% map[shared parameterID channelID]

lambda = parameters.lambda;
zemit0 = parameters.zemit0;
objStage0= parameters.objStage0; 
coder.varsize('zemit');
coder.varsize('zemitStack');
zemitStack = parameters.zemitStack';
zemit = zemitStack;
objStageStack = parameters.objStageStack;
sizeX = parameters.sizeX;
sizeY = parameters.sizeY; 
sizeZ = parameters.sizeZ;
XPupil=parameters.XPupil;
PolarizationVector=parameters.PolarizationVector;
Amplitude=parameters.Amplitude;
wavevector=parameters.wavevector;
wavevectorzmed=parameters.wavevectorzmed;
normfac=parameters.normfac;
allzernikes=parameters.allzernikes;
Ax=parameters.Ax;
Bx=parameters.Bx;
Dx=parameters.Dx;
Ay=parameters.Ay;
By=parameters.By;
Dy=parameters.Dy;
normIntensity = parameters.normIntensity;

numparams = parameters.numparams;
numAberrations = parameters.numAberrations;
ztype = parameters.ztype;

parameterID = numAberrations+5;
channelID = parameters.sizeZ;

repTheta = zeros(parameterID,channelID);

for i = 1:numparams
    if map(i,1)==0
        repTheta(map(i,2),map(i,3)) = theta(i);
    elseif map(i,1)==1||map(i,1)==2
        repTheta(map(i,2),1:channelID) = theta(i)*ones(1,channelID);
    end
end

% enlarging ROI in order to accommodate convolution with bead
% enlarging axial range in order to accommodate convolution with bead



if isfield(parameters,'bead')
  if parameters.bead == true
    ImageSizex = (parameters.pixelSizeX*parameters.sizeX/2);
    ImageSizey = (parameters.pixelSizeY*parameters.sizeY/2);
    DxImage = 2*ImageSizex/sizeX;
    DyImage = 2*ImageSizey/sizeY;
    beaddiameter = parameters.beaddiameter;
    DeltaMx = 2*ceil((beaddiameter/2)/DxImage);
    sizeX = sizeX + DeltaMx;
    ImageSizex = ImageSizex+DeltaMx*DxImage/2;
    DeltaMy = 2*ceil((beaddiameter/2)/DyImage);
    sizeY = sizeY + DeltaMy;
    ImageSizey = ImageSizey+DeltaMy*DyImage/2;
    
    NA = parameters.NA;
    Npupil = parameters.Npupil;
    ImageSizex = ImageSizex*NA/lambda;
    ImageSizey = ImageSizey*NA/lambda;
    PupilSize = 1.0;
    [Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,sizeX);
    [Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,sizeY);

    zemit = parameters.zemitStack;
    DzImage = zemit(1,end)-zemit(1,end-1);
    DeltaMz = 2*ceil((beaddiameter/2)/DzImage);
    sizeZ = sizeZ + DeltaMz;
    zmin = min(zemit)-DeltaMz*DzImage/2;
    zmax = max(zemit)+DeltaMz*DzImage/2;

    left_range = zmin:DzImage:min(zemit);
    right_range = max(zemit):DzImage:zmax;
    zemit = [left_range zemit(:,2:end-1) right_range];
    % zemit = zmin:DzImage:zmax;

  end
end


% calculation aberration function
Waberration = zeros(size(XPupil));
zernikecoefs = squeeze(repTheta(1:numAberrations,1));% seems to work only when aberration is shared??
zernikecoefs = normfac.*zernikecoefs;
for j = 1:numel(zernikecoefs)
  Waberration = Waberration+zernikecoefs(j)*squeeze(allzernikes(j,:,:));  
end

PhaseFactor = exp(2*pi*1i*Waberration/lambda);


numders = 3+numAberrations;

FieldMatrix = cell(2,3,sizeZ);
FieldMatrixDerivatives = cell(2,3,sizeZ,numders);

objStage = objStageStack(:)+objStage0(:);
if parameters.bead == true
  xemit = ones(1,sizeZ)* mean(repTheta(numAberrations+1,:));
  yemit = ones(1,sizeZ)* mean(repTheta(numAberrations+2,:));
  objStageStack = ones(1,sizeZ)* 0;
  if strcmp(ztype,'stage')
    objStage = zemit' + objStageStack(:)+objStage0(:);
    zemit = zemitStack(:)+zemit0(:);
  elseif strcmp(ztype,'emitter')
     objStage = objStageStack(:) + objStage0(:);
     zemit = (ones(1,sizeZ)*mean(repTheta(numAberrations+3,:))) + zemit+ zemit0(:);
  end
else
  xemit = repTheta(numAberrations+1,:);
  yemit = repTheta(numAberrations+2,:);
  if strcmp(ztype,'stage')
      objStage = repTheta(numAberrations+3,:)'+objStageStack(:)+objStage0(:);
      zemit = zemitStack(:)+zemit0(:);
  elseif strcmp(ztype,'emitter')
      objStage = objStageStack(:)+objStage0(:);
      zemit = mean(repTheta(numAberrations+3,:))+zemitStack(:)+zemit0(:);
    % zemit = zemitStack(:)+zemit0(:);
  end
end  



Nph = repTheta(numAberrations+4,:);
bg = repTheta(numAberrations+5,:);

% numparams = 26;
PSF = zeros(sizeX,sizeY,sizeZ);
PSFderivatives = zeros(sizeX,sizeY,sizeZ,numders);

for jz = 1:sizeZ
    for itel = 1:2
        for jtel = 1:3
            for jder = 1:numders
                FieldMatrixDerivatives{itel,jtel,jz,jder} = 1;
            end
        end
    end
end


for jz = 1:sizeZ
    % xyz induced phase contribution
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

                PupilFunction = 1i*wavevector{jder}.*PositionPhaseMask.*PupilMatrix;
                IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
                FieldMatrixDerivatives{itel,jtel,jz,numAberrations+jder} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

            end
            
            % pupil functions for z-derivative and FT to matrix elements

            if strcmp(ztype,'stage') % find stage position
                PupilFunction = 1i*wavevector{3}.*PositionPhaseMask.*PupilMatrix;
            elseif strcmp(ztype,'emitter') % find emitter position
                PupilFunction = 1i*wavevectorzmed.*PositionPhaseMask.*PupilMatrix;
            end
            IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
            FieldMatrixDerivatives{itel,jtel,jz,numAberrations+3} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

            
            % pupil functions for Zernike mode-derivative and FT to matrix elements
            %        n = 1;
            for jzer = 1:numAberrations
                
                PupilFunction = (2*pi*1i*normfac(jzer)*squeeze(allzernikes(jzer,:,:))/lambda).*PositionPhaseMask.*PupilMatrix;
                IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
                FieldMatrixDerivatives{itel,jtel,jz,jzer} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
                
            end
            
            
        end
    end
end


for jz = 1:sizeZ
    for jtel = 1:3
        for itel = 1:2
            PSF(:,:,jz) = PSF(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
            for jder = 1:numders

                PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) +...
                    (2/3)*real(conj(FieldMatrix{itel,jtel,jz}).*FieldMatrixDerivatives{itel,jtel,jz,jder});

            end
        end
    end
end

% 光子数归一化
PSF = PSF/normIntensity;
PSFderivatives = PSFderivatives/normIntensity;



% 3D convolution of the PSFs and derivatives with a bead
if isfield(parameters,'bead')
  if parameters.bead == true
    PSF = convn(bead,PSF,'same');
    [Mx,My,Mz] = size(PSF);
    tempderivs = zeros(Mx,My,Mz,numders);
    for jder = 1:numders
      tempderivs(:,:,:,jder) = convn(bead,squeeze(PSFderivatives(:,:,:,jder)),'same');
    end
    PSFderivatives = tempderivs;
  end
end



model = PSF;
sizeX = size(PSF,1);
sizeY = size(PSF,2);
sizeZ = size(PSF,3);
for i = 1:sizeZ
    model(:,:,i) = Nph(i)*PSF(:,:,i)+bg(i);
end

dudt = zeros(sizeX,sizeY,sizeZ,numAberrations+5);
for i = 1:sizeZ
    for jp = 1:2
        if shared(numAberrations+jp)~=2
            dudt(:,:,i,numAberrations+jp) = -1*Nph(i)*PSFderivatives(:,:,i,numAberrations+jp);
        end
    end
    for jp = 1:numAberrations
        if shared(jp)~=2
            dudt(:,:,i,jp) = Nph(i)*PSFderivatives(:,:,i,jp);
        end  
    end
    if shared(numAberrations+3)~=2
        dudt(:,:,i,numAberrations+3) = Nph(i)*PSFderivatives(:,:,i,numAberrations+3);
    end
end


dudt(:,:,:,numAberrations+4) = PSF;
dudt(:,:,:,numAberrations+5) = ones(size(PSF));




end