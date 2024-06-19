% pay attention to modify the initial parameters
parameters = set_parameters(Npixels, numel(zemit), zemit, standard);   % % The result of VISPR is affected by initial parameters===================================================================================================
parameters.aberrations(1,3) = -60;  % Vertical Astig
parameters.aberrations(2,3) = 0; % Oblique Astig
parameters.bead = true;
if exist( 'parameters_real', 'var') == 1 % determine simulation or experiment 
   parameters.NA = parameters_real.NA;
   parameters.refmed = parameters_real.refmed;
   parameters.refcov = parameters_real.refcov;
   parameters.refimm = parameters_real.refimm;
   parameters.lambda = parameters_real.lambda;
   parameters.pixelSizeX = parameters_real.pixelSizeX;
   parameters.pixelSizeY = parameters_real.pixelSizeY;
end

ref_PSF = gen_vectorial_psf(parameters);
% imageslicer(single(permute(ref_PSF,[2 1 3])));
% title('first reference')

P = [];  

for iter = 1 : iteration

    disp(['Iteration: ' num2str(iter) ' ...']);
    
    % Generate reference Z-postions' PSFs from pupil function
    if isempty(P) == 0
        parameters.aberrations(:,3) = P(1:zernikenum);
        ref_PSF = gen_vectorial_psf(parameters);
        ref_PSF = ref_PSF.* Nph + bg;
    else
        ref_PSF = ref_PSF.* Nph + bg;
    end
    
    % classify single molecules to giving reference Z-position images, and 
    %  average these Z-postion images in a certain group 
    %  calculate similartiy between single molecules and reference images, do
    %  X-Y registration
    
    empupil.imsz = size(ref_PSF,1);  % image的pixel大小
    [ims_Zcal_ave_plane1, index_record_Zplanes] = classify_onePlane_par(subregion_ch1, ref_PSF, empupil);
    
    % Use phase retrival to refine the pupil function and estimate the aberration using average 
    %  images in different Z postions
    disp(['Z-position index: ' num2str(index_record_Zplanes')]);
    now_zemit = zemit(index_record_Zplanes==1); 
    ims_Zcal_ave_plane1_sel = ims_Zcal_ave_plane1(:,:,index_record_Zplanes==1);

    if parameters.bead == true && mod(size(now_zemit,2),2)==0  % Regarding the circumstance of beads
        now_zemit(:,end) = [];
        ims_Zcal_ave_plane1_sel(:,:,end) = [];
    end
    now_zmin = min(now_zemit);
    now_zmax = max(now_zemit);
    now_zrange = now_zmax - now_zmin;
    if now_zrange == 0   
        error('Error: please check out the setting of initial parameters');  
    end
    disp(['zmin: ',num2str(now_zmin),';  ','zmax: ',num2str(now_zmax),';  ','zrange: ', num2str(now_zrange)])

    paraFit = set_parameters_fit(ims_Zcal_ave_plane1_sel, now_zemit, standard);  % pay attention to modify the initial parameters such as refrective index, NA===================================================================================================
    paraFit.bead = true;
    if exist( 'parameters_real', 'var') == 1 
       paraFit.NA = parameters_real.NA;
       paraFit.refmed = parameters_real.refmed;
       paraFit.refcov = parameters_real.refcov;
       paraFit.refimm = parameters_real.refimm;
       paraFit.lambda = parameters_real.lambda;
       paraFit.pixelSizeX = parameters_real.pixelSizeX;
       paraFit.pixelSizeY = parameters_real.pixelSizeY;
    end

    [P,model,err] = mine_MLE_fit(ims_Zcal_ave_plane1_sel, paraFit.thetainit, paraFit, paraFit.shared, 0.1,zernikenum);
    

    disp('...PR Successfully finished!');
    display(['Updated Zernike(1-16): ', num2str(P(1:16)')]);
    if iter~=iteration
        close all
    end

end



%% plot results
if exist('Zernike','var') == 0
    f1 = figure;
    axMode =axes(f1);
    bar(axMode, P(1:zernikenum),0.3);
    legend(axMode,'Fitting Result')
    orders = paraFit.aberrations(:,1:2);
    for k=size(orders):-1:1
        axn{k}=[num2str(orders(k,1)) ',' num2str(orders(k,2))];
    end
    axMode.XTick=1:length(axn);
    axMode.XTickLabel=axn;
    ylim([-70 50]);
    ylabel('Zernike rms value(nm)') ,xlabel('Zernike mode')
else
    f1 = figure;
    aberrations = Zernike;
    axMode =axes(f1);
    bar(axMode,[aberrations(:,3), P(1:zernikenum)]);
    legend(axMode,{'Ground truth','Vectorial INSPR'})
    orders = paraFit.aberrations(:,1:2);
    for k=size(orders):-1:1
        axn{k}=[num2str(orders(k,1)) ',' num2str(orders(k,2))];
    end
    axMode.XTick=1:length(axn);
    axMode.XTickLabel=axn;
    xtickangle(45)
    ylabel('Zernike rms value(nm)','FontWeight','bold') ,xlabel('Zernike mode','FontWeight','bold')
    disp('...Successfully finished!');
end



if exist('Zernike','var') == 0
    PupilSize = 1.0;
    Npupil = paraFit.Npupil;
    DxyPupil = 2*PupilSize/Npupil;
    XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
    [YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
    ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
    Waberration_fit = zeros(size(XPupil));
    orders = paraFit.aberrations(:,1:2);
    zernikecoefs_fit = P(1:zernikenum);
    normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
    zernikecoefs_fit = normfac.*zernikecoefs_fit;
    allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
    for j = 1:numel(zernikecoefs_fit)
        Waberration_fit = Waberration_fit+zernikecoefs_fit(j)*squeeze(allzernikes(j,:,:));
    end
    Waberration_fit = Waberration_fit.*ApertureMask; % nm
    f2 = figure;
    axpupil = axes(f2);
    imagesc(Waberration_fit);title('Fit pupil');
    ylabel('Y(pixel)') ,xlabel('X(pixel)')
    c_title=colorbar;
%     caxis([-250 600]);
    title(c_title,'nm')

else
    PupilSize = 1.0;
    Npupil = paraFit.Npupil;
    DxyPupil = 2*PupilSize/Npupil;
    XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
    [YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
    ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
    
    Waberration_fit = zeros(size(XPupil));
    Waberration_true = zeros(size(XPupil));
    orders = paraFit.aberrations(:,1:2);
    zernikecoefs_fit = P(1:zernikenum);
    zernikecoefs_true = aberrations(:,3);
    
    normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
    zernikecoefs_fit = normfac.*zernikecoefs_fit;
    zernikecoefs_true = normfac.*zernikecoefs_true;
    allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
    for j = 1:numel(zernikecoefs_fit)
      Waberration_fit = Waberration_fit+zernikecoefs_fit(j)*squeeze(allzernikes(j,:,:));
      Waberration_true = Waberration_true+zernikecoefs_true(j)*squeeze(allzernikes(j,:,:));
    end
    
    Waberration_fit = Waberration_fit.*ApertureMask; % nm
    Waberration_true = Waberration_true.*ApertureMask; % nm
    % pupil = Waberration_fit/paraFit.lambda*2*pi;  % rad
    RMSE_value = RMSE(Waberration_fit,Waberration_true);
    disp(['RMSE_value = ',num2str(RMSE_value),' nm'])
    Waberration_diff = Waberration_true - Waberration_fit;
    
    ApertureMask(ApertureMask==0) = NaN;
    Waberration_fit = Waberration_fit.*ApertureMask;
    Waberration_true = Waberration_true.*ApertureMask;
    Waberration_diff = Waberration_diff.*ApertureMask;
    
    f2 = figure;
    set(gcf,'unit','normalized','position',[0.1,0.1,0.75,0.33]);
    axpupil = axes(f2);
    min_ab = min(Waberration_true(:));
    max_ab = max(Waberration_true(:));
    
    % imagesc(axpupil,pupil);
    subplot(1,3,1);
    h = imagesc(Waberration_fit);
    set(h,'alphadata',~isnan(Waberration_fit));
    % axis(axpupil,'equal')
    % axis(axpupil,'tight')
    % ylabel('Y(pixel)') ,xlabel('X(pixel)')
    axis off
    title('Vectorial insitu Pupil','FontWeight','bold');
    % c_title=colorbar;
    % caxis([min_ab max_ab]);
    % title(c_title,'nm')
    
    subplot(1,3,2);
    h = imagesc(Waberration_true);
    set(h,'alphadata',~isnan(Waberration_true));
    % axis(axpupil,'equal')
    % axis(axpupil,'tight')
    title('Ground truth','FontWeight','bold');
    % ylabel('Y(pixel)') ,xlabel('X(pixel)')
    axis off
    % c_title=colorbar;
    % title(c_title,'nm')

    subplot(1,3,3);
    h = imagesc(Waberration_diff);
    set(h,'alphadata',~isnan(Waberration_diff));
    % axis(axpupil,'equal')
    % axis(axpupil,'tight')
    title('Residual error','FontWeight','bold');
    % ylabel('Y(pixel)') ,xlabel('X(pixel)')
    axis off
    
    c=colorbar;
    % title(c,'nm')
    colormap('jet');
    caxis([min_ab max_ab]);
    
    set(c,'YTick',[min_ab max_ab]);
    min_lambda = roundn(min_ab/parameters_real.lambda,-1);
    max_lambda = roundn(max_ab/parameters_real.lambda,-1);
    set(c,'YTickLabel',{[num2str(min_lambda),'\lambda'],[num2str(max_lambda),'\lambda']}) 
    set(c,'FontName','time','FontSize',14,'FontWeight','bold'); 
%     title(c_title,'nm')
    
    

end