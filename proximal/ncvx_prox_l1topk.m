%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_l1topk
clc;clear all;close all;
% min_t 0.5*alpha*t^2 + beta*t + lambda ||x+tei||_1 - lambda || x+tei ||_{top-k}

for iter = 1:100000000
    
    rand('seed',iter);
    randn('seed',iter);
    m = 10;
    
    alpha = rand(1)*100*rand(1);
    beta = randn(1)*100*rand(1);    beta = beta*max(0,randn(1));
    x = randn(m,1)*100*rand(1);    x = x.*max(0,randn(m,1));
    lambda = rand(1)*100*rand(1);     lambda = lambda*max(0,randn(1));
    i = randperm(m,1);
    k = randperm(m,1);
    
    HandleObj = @(t)ComputeObj(t,alpha,beta,lambda,x,i,k);
    x1 = fminsearch(HandleObj,0);
    x2 = ncvx_prox_l1topk_(alpha,beta,lambda,x,i,k);
    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    fprintf('iter:%.6d, fobj:%.1e %.1e, diff: %.5e\n',iter,f1,f2,f1-f2);
    if(f2>f1+1e-5*abs(mean([f1;f2])))
        f1
        f2
        x1
        x2
        dddd
    end
    
    
end
%}



function best_t = ncvx_prox_l1topk_(alpha,beta,lambda,x,i,k)
% min_t 0.5 alpha t^2 + beta t + lambda ||x+tei||_1 - lambda || x+tei ||_{top-k}
m = length(x);ei = zeros(m,1); ei(i)=1;
ts = [-x(i);-beta/alpha;(lambda-beta)/alpha;(-beta-lambda)/alpha];
HandleObj = @(t)0.5*alpha*t^2 + beta*t + lambda*abs(x(i)+t) - lambda*topksum(x+t*ei,k);
best_fobj = inf;
for i = 1:length(ts)
    test_obj = HandleObj(ts(i));
    if(test_obj<best_fobj)
        best_fobj = test_obj;
        best_t = ts(i);
    end
end

function r = topksum(x,k)
val = sort(abs(x),'descend');
r = sum(val(1:k));


function [f] = ComputeObj(t,alpha,beta,lambda,x,i,k)
x1 = x; x1(i) = x(i)+t;
f = 0.5*alpha*t*t + beta*t + lambda*norm(x1,1) - lambda*topksum(x1,k);
