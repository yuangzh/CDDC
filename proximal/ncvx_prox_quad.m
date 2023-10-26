%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_quad_bound
clc;clear all;close all;
% 0.5*a*t*t + b*t, s.t. l <= t<= u

for iter = 1:100000000
    m = 1;
    a = max(0,randn(m,1))*100*randn(1);
    b = max(0,randn(m,1))*100*randn(1);
    lu = sort(randn(2,1)*100.*max(0,rand(2,1)));
    l = lu(1);
    u = lu(2);
    
    proj = @(t) max(min(t,u),l);
    HandleObj = @(t) 0.5*a*proj(t)^2 + b*proj(t);
    
    x1 = fminbnd(HandleObj,l,u);
    x2 = ncvx_prox_quad___(a,b,l,u);
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





function [ts] = ncvx_prox_quad___(alpha,beta,l,u)
% 0.5*alpha*t*t + beta*t, s.t. l <= t<= u
% Note that {alpha,beta,l,u} can be vectors

handle_f = @(t)0.5.*alpha.*t.*t + beta.*t;
proj = @(x)max(min(x,u),l);
x1 = l;
x2 = u;
x3 = proj(-beta./alpha);
% x_i = arg min_{t \in {x1(i),x2(i),x3(i)} }  f(t)
% with i = 1,..., n
% handle_f: R => R
% x1: n x 1
% x2: n x 1
% x3: n x 1
dim = length(x1);
num = 3;
X = [x1 x2 x3]';
f1 = handle_f(x1);
f2 = handle_f(x2);
f3 = handle_f(x3);
F = [f1 f2 f3];
[~,I] = min(F,[],2);
x = X(I + [1:num:dim*num]'-1);
ts = reshape(x,dim,1);


function [t] = ncvx_prox_quad_bound_scalar(alpha,beta,l,u)
% 0.5*alpha*t*t + beta*t, s.t. l <= t<= u
f = @(t)0.5*alpha*t*t + beta*t;
proj = @(x)max(min(x,u),l);
t1 = l; t2 = u; t3 = proj(-beta./alpha);

f1 = f(t1); f2 = f(t2); f3 = f(t3);
if(f1<=f2 && f1<= f3)
    t = t1;
elseif(f2<=f1 && f2<= f3)
    t = t2;
else
    t = t3;
end
