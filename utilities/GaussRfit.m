function [I,sigma,bg,fit_im]=GaussRfit(obj,startpoint,input_im)

R1=obj.OTFratioSize;
R=obj.PSFsize;
estimate=fminsearch(@(x) Gauss2(x,input_im,R1,R,obj.Pixelsize),startpoint,optimset('MaxIter',50,'Display','off'));

I=estimate(:,1);
tmp=estimate(:,2);
tmp(tmp>5)=5;
sigma = tmp;

bg=0;
scale=R*obj.Pixelsize;
[xx,yy]=meshgrid(-R/2:R/2-1,-R/2:R/2-1);

X=abs(xx)./scale;
Y=abs(yy)./scale;
% fit_im=I.*exp(-X.^2./2./sigmax^2).*exp(-Y.^2./2./sigmay^2)+bg;
fit_im=I.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2)+bg;
end