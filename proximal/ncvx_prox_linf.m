%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_linf
clc;clear all;close all;
% 0.5 alpha t*t + beta*t  - lambda || t a + b||_inf

for iter = 1:100000000
    n = 3;
    lambda = rand(1)*100*rand(1);
    alpha = rand(1)*100*rand(1);
    beta = randn(1)*100*rand(1);
    a = randn(n,1)*100*rand(1).*max(0,randn(n,1));
    b = randn(n,1)*100*rand(1).*max(0,randn(n,1));
    HandleObj = @(t) 0.5*alpha*t*t + beta*t - lambda*max(abs(t*a+b));
    x1 = fminsearch(HandleObj,0);
    x2 = ncvx_prox_linf___(alpha,beta,lambda,a,b);
    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    fprintf('iter:%.6d, fobj:%.1e %.1e, diff: %.5e\n',iter,f1,f2,f1-f2);
    if(f2>f1+1e-5*abs(mean([f1;f2])))
        x1
        x2
        dddd
    end
end

%}

function [t] = ncvx_prox_linf___(alpha,beta,lambda,a,b)
% 0.5 alpha t*t + beta*t  - lambda || t a + b||_inf
% Computation Complexity: O(n)

ts = [-beta/alpha + lambda/alpha*a;  -beta/alpha - lambda/alpha*a];
aa = [a;a];
bb = [b;b];
for i = 1:length(ts)
    t  = ts(i);
    his(i) = 0.5*alpha*t*t + t*beta - lambda*(abs(t*aa(i)+bb(i)));
end
[val,ind] = min(his);
t = ts(ind);

