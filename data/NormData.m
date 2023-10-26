
function [x] = NormData(x,type)
% x : N x d
%  normalize each column of x
% type
% 1: normalize each column of x to have unit norm
% 2: normalize each column of x to have zero mean and one standard deviation
% 3: Scale each column of x to [-1 1]
if(type==1)
    d = size(x,2);
    feaNorm = max(1e-9,full(sum(x.^2,1))');
    x = x*spdiags(feaNorm.^-.5,0,d,d);
elseif(type==2)
    d=size(x,2);
    tol=1e-5;
    meanxapp=zeros(1,d); stdxapp=zeros(1,d);
    for i=1:d,
        meanxapp(i)=mean(x(:,i));
        stdxapp(i)=std(x(:,i));
        if(stdxapp(i)<tol),stdxapp(i)=1;end;
        x(:,i)= (x(:,i) - meanxapp(i)) / stdxapp(i) ;
    end;
%     x=zscore(x);
elseif(type==3)
    low = -1;up = 1;
    len = up - low;
    cmin = min(x,[],1);
    cmax = max(x,[],1);
    rowsNum = size(x,1);
    maxSubMin = cmax - cmin;
    DataMin = repmat(cmin,rowsNum,1);
    DmaxSubMin = repmat(maxSubMin,rowsNum,1);
    x = low+(x-DataMin)./DmaxSubMin * len;
end
 


