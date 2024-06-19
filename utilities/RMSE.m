function varargout= RMSE(varargin)

if nargin<2
    error('输入参数必须大于等于2.');
elseif nargin==2
    varargout{1}=h(varargin{1},varargin{2});
end
end


function r=h(f,g)
%函数功能：
%      函数r=h(f,g)求出两幅图像的均方根误差r
%输入参数：
%      f----标准图像
%      g----融合后的图像
%-------------------------------------%

f=double(f);
g=double(g);
[m,n]=size(f);

temp=[];
for i=1:m
    for j=1:n
        temp(i,j)=(f(i,j)-g(i,j))^2;
    end
end

r=sqrt(sum(sum(temp))/(m*n));

end
