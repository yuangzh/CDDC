%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_linfl2
clc;clear all;close all;
% min_t 0.5*alpha*t^2 + beta t - lambda || x+tei ||, s.t. lb <= t <= ub

for iter = 1:100000000
    n = 3;
    lambda = rand(1)*100*rand(1);
    alpha = rand(1)*100*rand(1);
    beta = randn(1)*100*rand(1);
    x = randn(n,1)*100*rand(1).*max(0,randn(n,1)); 
    lbub = randn(2,1)*100*rand(1).*max(0,randn(2,1)); lbub = sort(lbub);
    lb = lbub(1);
    ub = lbub(1);
    i = randperm(n,1);
    ei = zeros(n,1); ei(i)=1;
    proj = @(t) max(min(t,ub),lb);

    HandleObj = @(t) 0.5*alpha*proj(t)^2 + beta*proj(t) - lambda*norm(x+proj(t)*ei);

    x1 = fminbnd(HandleObj,-10000,10000);
    x2 = ncvx_prox_linfl2___(alpha,beta,lambda,x,i,lb,ub);
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

function best_t = ncvx_prox_linfl2___(alpha,beta,lambda,x,i,lb,ub)
% min_t 0.5*alpha*t^2 + beta t - lambda || x+tei ||, s.t. lb <= t <= ub
tt =  ncvx_prox_l2(alpha,beta,lambda,x,i);
tt = max(lb,min(tt,ub));
ts = [lb;ub;tt];
HandleObj = @(t)ComputeObj(t,alpha,beta,lambda,x,i,lb,ub);
for i=1:length(ts)
    fs(i) = HandleObj(ts(i));
end
[~,j] = min(fs);
best_t = ts(j);



function [ts] = ncvx_prox_l2(alpha,beta,lambda,x,i)
% Compute all critical points
% 0.5 alpha t^2 + beta*t - lambda sqrt ( ||  x+tei ||_2^2 )
% 0.5 alpha t^2 + beta*t - lambda sqrt ( ||  x+tei ||_2^2 )
% 0.5 alpha t^2 + beta*t - sqrt ( a0 + a1*t + a2*t*t )
% grad = alpha*t + beta - (a0 + a1*t + a2*t*t)^(-1/2)*(0.5*a1+a2*t) = 0
% alpha*t + beta = (a0 + a1*t + a2*t*t)^(-1/2)*(0.5*a1+a2*t)
% (alpha*t + beta)^2 * (a0 + a1*t + a2*t*t)  = (0.5*a1+a2*t)^2

a0 = lambda*lambda*x'*x;
a1 = lambda*lambda*2*x(i);
a2 = lambda*lambda;

% t = 1.234;
% r1 = (alpha*t + beta)^2 * (a0 + a1*t + a2*t*t)  - (0.5*a1+a2*t)^2
c4 = alpha*alpha*a2;
c3 = alpha*alpha*a1 + 2*alpha*beta*a2 ;
c2 = alpha*alpha*a0   +  2*alpha*beta*a1  + beta*beta*a2 - a2*a2;
c1 = 2*alpha*beta*a0 + beta*beta*a1 - 0.5*a1*a2 - a2*0.5*a1;
c0 = beta*beta*a0  - 0.5*a1*0.5*a1;

% r2 = @(t) c4*t^4 + c3*t^3 + c2*t^2 + c1*t + c0
ts = roots([c4;c3;c2;c1;c0]);
ts = [real(ts);-beta/alpha];

function [f] = ComputeObj(t,alpha,beta,lambda,x,i,lb,ub)
% min_t 0.5*alpha*t^2 + beta*t - lambda || x+tei ||_2
% s.t. lb <= t <= ub
t = max(lb,min(t,ub));
xx = x'*x;
f = 0.5*alpha*t^2 + beta*t - lambda*sqrt( xx + 2*x(i)*t + t*t );