%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_l1_slow
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
    x2 = ncvx_prox_l1_slow___(alpha,beta,lambda,a,b);
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

function best_t = ncvx_prox_l1_slow___(alpha,beta,lambda,a,b)
% min_t 0.5*alpha*t^2+beta*t-lambda*||ta+b||_1
I=find(a==0); a(I)=[]; b(I)=[];
J=find(a<0); a(J)=-a(J); b(J)=-b(J);
e=-b./a; [~,o]=sort(e,'ascend');
inf=1e200; z=[-inf;e(o);inf]; his=[]; ts=[];
for i = 1:(length(z)-1)
  tt=(z(i)+z(i+1))/2; 
  t=-beta/alpha+lambda/alpha*sign(tt*a+b)'*a;
  if((t>=z(i)-1e-16)&&(t<=z(i+1)+1e-16))
     ft=0.5*alpha*t*t+beta*t-lambda*norm(t*a+b,1); 
     his=[his;ft]; ts=[ts;t];   end
end
[~,ind]=min(his); best_t = ts(ind);
