
function [sse,Model]=Gauss2(x,input_im,R1,R,pixelsize)
I=x(1);
% sigmax=x(2);
% sigmay=x(3);
sigma=x(2);
bg=0;
[xx,yy]=meshgrid(-R1/2:R1/2-1,-R1/2:R1/2-1);

scale=R*pixelsize;
X=abs(xx)./scale;
Y=abs(yy)./scale;
% Model=I.*exp(-X.^2./2./sigmax^2).*exp(-Y.^2./2./sigmay^2)+bg;
Model=I.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2)+bg;

sse=sum(sum((Model-input_im).^2));
end