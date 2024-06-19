%%  1. Set parameters
rep = 100;  % Number of repetitions per position
x_pixelsize = parameters_real.pixelSizeX;   
z_pixelsize = 10;
Nphoton = 5000;    
bg = 30;
roi = 25;
r_roi = (roi-1)/2;

zrange = -600:20:600;  % nm
znum = size(zrange,2);
x0 = 0.000000;
y0 = 0.000000;

mid = r_roi+1;
raw_result = zeros(rep, znum, 3);  

std_compare  = zeros(znum, 6);  % xyz localization precision
accuracy_compare  = zeros(znum, 6);  % xyz localization accuracy
mean_compare = zeros(znum, 6);
crlb = zeros(znum, 3); 


%% 2. Calculate the localization accuracy and precision for each depth
tic

data_cell = cell(znum,1);

for index = 1:znum   
    raw_data = zeros(roi, roi, rep); 
    zz = size(coeff,3)/2;
    [x, y, z, ~] = size(coeff);
    off = (x - roi-1)/2;
    z_set = zrange(index)/z_pixelsize+zz;

    xc = x0;
    yc = y0;
    zc = z_set;

    delta_x = -1*xc;
    xstart = floor(delta_x);
    delta_x = delta_x - xstart;

    delta_y = -1*yc;
    ystart = floor(delta_y);
    delta_y = delta_y - ystart;

    delta_z = zc - floor(zc);
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2((delta_x),(delta_y),(delta_z));

    z_index = floor(zc);


    if z_index <=0
        z_index = 1;
    end
    if z_index >z
        z_index = z;
    end

    test1 = zeros(roi, roi);
    for ii = 0:roi-1
        for jj = 0:roi-1
            temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,z,delta_f,coeff);
            test1(ii+1,jj+1)=temp;
        end
    end
    total_N = sum(test1(:));
    N =  Nphoton/total_N;
    total_N = N*sum(test1(:));
    parfor i = 1:rep
        raw_data(:, :, i) = N* test1;
        
    end

  
    ims = zeros(size(raw_data));

%     tic
    parfor i = 1:rep
        tmp = raw_data(:,:,i);
        tmp = uint16(tmp)+bg;
        ims(:,:,i) = double(imnoise(tmp, 'poisson'));
    end
%     toc
    data_cell{index,1} = ims;   
end

Real_mean_x = zeros(znum, 1);
Real_mean_y = zeros(znum, 1);
Real_mean_z = zeros(znum, 1);

VISPR_mean_x = zeros(znum, 1);
VISPR_mean_y = zeros(znum, 1);
VISPR_mean_z = zeros(znum, 1);



for index = 1:znum   
    ims = data_cell{index,1};
    sCMOSvarmap = 0;
    [result,CRLB] = mleFit_LM( single(ims), 6, 50, single(coeff_real), sCMOSvarmap,1);
    % result: X, Y, Photons, Background, Z, Iterations
    x_real = result(:,1);
    y_real = result(:,2);
    z_real = result(:,end-1);
    crlb_x = mean(CRLB(:,1))*x_pixelsize;
    crlb_y = mean(CRLB(:,2))*x_pixelsize;
    crlb_z = mean(CRLB(:,end))*z_pixelsize;
    
    [result_VISPR,CRLB_XX] = mleFit_LM( single(ims), 6, 50, single(coeff_VISPR), sCMOSvarmap,1);
    x_VISPR = result_VISPR(:,1);
    y_VISPR = result_VISPR(:,2);
    z_VISPR = result_VISPR(:,end-1);

    Real_mean_x(index,1) = mean(x_real);
    Real_mean_y(index,1) = mean(y_real);
    Real_mean_z(index,1) = mean(z_real);

    VISPR_mean_x(index,1) = mean(x_VISPR);
    VISPR_mean_y(index,1) = mean(y_VISPR);
    VISPR_mean_z(index,1) = mean(z_VISPR);


end

% Calculate global bias
% myfun = @(x,xdata) x + xdata;
% delta_x_VISPR = lsqcurvefit(myfun,0,VISPR_mean_x, Real_mean_x);
% delta_y_VISPR = lsqcurvefit(myfun,0,VISPR_mean_y, Real_mean_y);
% delta_z_VISPR = lsqcurvefit(myfun,0,VISPR_mean_z,Real_mean_z);

for index = 1:znum   
    ims = data_cell{index,1};
    sCMOSvarmap = 0;
    z_set = zrange(index)/z_pixelsize+zz;

    [result,CRLB] = mleFit_LM( single(ims), 6, 50, single(coeff_real), sCMOSvarmap,1);
    x_real = result(:,1);
    y_real = result(:,2);
    z_real = result(:,end-1);
    crlb_x = mean(CRLB(:,1))*x_pixelsize;
    crlb_y = mean(CRLB(:,2))*x_pixelsize;
    crlb_z = mean(CRLB(:,end))*z_pixelsize;
    
    [result_VISPR,CRLB_XX] = mleFit_LM( single(ims), 6, 50, single(coeff_VISPR), sCMOSvarmap,1);
    x_VISPR = result_VISPR(:,1);
    y_VISPR = result_VISPR(:,2);
    z_VISPR = result_VISPR(:,end-1);

    raw_result(:,index,1) = x_VISPR;
    raw_result(:,index,2) = y_VISPR;
    raw_result(:,index,3) = z_VISPR;

    std_compare(index,1) = std(x_real)*x_pixelsize;
    std_compare(index,2) = std(y_real)*x_pixelsize;
    std_compare(index,3) = std(z_real)*z_pixelsize;
    std_compare(index,4) = std(x_VISPR)*x_pixelsize;
    std_compare(index,5) = std(y_VISPR)*x_pixelsize;
    std_compare(index,6) = std(z_VISPR)*z_pixelsize;
    
    accuracy_compare(index,1) = sqrt(mean((x_real - (x0+mid)).^2))*x_pixelsize;
    accuracy_compare(index,2) = sqrt(mean((y_real - (y0+mid)).^2))*x_pixelsize;
    accuracy_compare(index,3) = sqrt(mean((z_real - z_set).^2))*z_pixelsize;
    accuracy_compare(index,4) = sqrt(mean((x_VISPR - (x0+mid)).^2))*x_pixelsize;
    accuracy_compare(index,5) = sqrt(mean((y_VISPR - (y0+mid)).^2))*x_pixelsize;
    accuracy_compare(index,6) = sqrt(mean((z_VISPR - z_set).^2))*z_pixelsize;
    
    crlb(index,1) = crlb_x;
    crlb(index,2) = crlb_y;
    crlb(index,3) = crlb_z;

    mean_compare(index,1) = (mean(result(:,1))-mid)*x_pixelsize;
    mean_compare(index,2) = (mean(result(:,2))-mid)*x_pixelsize;
    mean_compare(index,3) = (mean(result(:,end-1))-zz)*z_pixelsize;
    mean_compare(index,4) = (mean(result_VISPR(:,1))-mid)*x_pixelsize;
    mean_compare(index,5) = (mean(result_VISPR(:,2))-mid)*x_pixelsize;
    mean_compare(index,6) = (mean(result_VISPR(:,end-1))-zz)*z_pixelsize;

end
toc




%% 3. Show results
% close all
pointsz = 14;
[all_themes, all_colors] = GetColors();
% % Show precision
figure;
subplot(331)
plot(zrange, std_compare(:,1), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, std_compare(:,4), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization precision in x (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')



subplot(332)
plot(zrange, std_compare(:,2), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, std_compare(:,5), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization precision in y (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')



subplot(333)
plot(zrange, std_compare(:,3), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, std_compare(:,6), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization precision in z (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')


% Show accuracy
% figure
subplot(334)
plot(zrange, accuracy_compare(:,1), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, accuracy_compare(:,4), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization accuracy in x (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')


% figure
subplot(335)
plot(zrange, accuracy_compare(:,2), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, accuracy_compare(:,5), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization accuracy in y (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')


%figure
subplot(336)
plot(zrange, accuracy_compare(:,3), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, accuracy_compare(:,6), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF ','FontSize',12,'location','north');
ylabel('localization accuracy in z (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')



% Show mean value
%figure;
subplot(337)
plot(zrange, mean_compare(:,1), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, mean_compare(:,4), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization mean value in x (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')


% figure
subplot(338)
plot(zrange, mean_compare(:,2), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, mean_compare(:,5), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization mean value in y (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')


%figure
subplot(339)
plot(zrange, mean_compare(:,3), 'o-', 'LineWidth', 2,'MarkerSize', 3);hold on
plot(zrange, mean_compare(:,6), '*--', 'LineWidth', 2,'MarkerSize', 3);
legend('Real PSF', 'VISPR PSF','FontSize',12,'location','north');
ylabel('localization mean value in z (nm)','FontSize',11,'FontWeight','bold')
xlabel('z (nm)','FontSize',11,'FontWeight','bold')





