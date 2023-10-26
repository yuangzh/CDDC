%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_l
clc;clear all;close all;
rand('seed',0);
randn('seed',0);
% min_t 0.5*alpha*t^2+beta*t-lambda*||tg+d||_1

for iter = 4:100000000
    m = 10;
    beta = randn(1)*100*rand(1);
    a = randn(m,1)*100*rand(1).*max(0,randn(m,1));
    b = randn(m,1)*100*rand(1).*max(0,randn(m,1));
    lambda = rand(1)*100*rand(1);
    alpha = rand(1)*100*rand(1);
    
    HandleObj = @(t) 0.5*alpha*t^2 + beta*t - lambda*norm(t*a+b,1);
    x1 = fminsearch(HandleObj,0);
    x2 = ncvx_prox_l1___(alpha,beta,lambda,a,b);
    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    fprintf('iter:%.6d, fobj:%.1e %.1e, diff: %.5e\n',iter,f1,f2,f1-f2);
    if(f2>f1+1e-5*abs(mean([f1;f2])))
        f1
        f2
        x1
        x2
        error('too large!')
    end
end

%}



function best_t = ncvx_prox_l1___(alpha,beta,lambda,a,b)
% min_t 0.5*alpha*t^2+beta*t-lambda*||ta+b||_1
I=find(a==0);a(I)=[];b(I)=[];
J=find(a<0);a(J)=-a(J); b(J)=-b(J);
e = -b./a; [~,o] = sort(e,'ascend'); a = a(o); b = b(o);
hf = @(t)0.5*alpha*t^2+t*beta-lambda*norm(t*a+b,1);
z = [-1e200;e(o);1e200]; ts = -beta/alpha; fobjs = hf(ts); 
for i = 1:(length(z)-1)
    if(i==1),r=-sum(a);
    else r=r+2*a(i-1);   end
    tt=(lambda*r-beta)/alpha;
    if((tt>=z(i)-1e-16)&&(tt<=z(i+1)+1e-16))
        fobjs = [fobjs;hf(tt)]; ts = [ts;tt];   end
end
[~,j]=min(fobjs); best_t = ts(j);

