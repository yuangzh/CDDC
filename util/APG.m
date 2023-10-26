%{
function test_APG
clear, clc;
%  min  1/2 || A x - y||^2 + lambda * ||x||_1
% f(x) + g(x)
% 0.5 L ||x-xt||^2 + <x-xt,g> + g(x)
% 0.5 L ||x-(xt-g/L)||^2 + g(x)
% proximal mapping:

% 0.5 theta ||x - a||^2 + g(x)


m=1000;  n=100;    % The data matrix is of size m x n
A=randn(m,n);       % the data matrix
y = randn(m,1);
lambda=0.2;
HandleObjSmooth = @(x)computeObj(x,A,y);
HandleObjNonSmooth = @(x)lambda*sum(abs(x));
x=zeros(n,1);
HandleProx = @(a,theta)computeprox(a,theta,lambda);
[x1, his]= APG___(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,100);
plot(his)

function [fobj,grad] = computeObj(x,A,y)
diff = A*x-y;
fobj = 1/2*norm(diff)^2 ;
grad = A'*diff ;

function [x] = computeprox(a,theta,lambda)
% 0.5 theta ||x - a||^2 + g(x)
[x] = threadholding_l1(a,lambda/theta);

function [x] = threadholding_l1(a,lambda)
% solving the following OP:
% min_{x} 0.5 ||x - a||^2 + lambda * sum(abs(x))
x = sign(a).*max(0,abs(a)-lambda);
%}


function [x_best,his] = APG___(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,maxiter)
% This program solves the following optimization problem:
% f(x) + g(x)
% where we assume that f is smooth g is non-smooth
% HandleObjSmooth:           x   ->  [fobj,grad]
% HandleObjNonSmooth:        x   ->  [fobj]
% HandleProx:          [a,theta] ->  arg min_{x} 0.5 theta || x - a ||^2 + g(x)


L = 0.1; % initial lipschitz parameter
[n,d]=size(x);
flag=0;
xp = x; xxp=zeros(n,d);
alpha=1; s=x; x_best = x;
f_best = HandleObjSmooth(x)+ HandleObjNonSmooth(x); his = [f_best];
for iter=1:maxiter,
    [f_old,g_old] = HandleObjSmooth(s);
    xp=x;
    for line_search =1:10000,
        v=s-g_old/L;
        % min_{x} 0.5 L ||x-(xt-g/L)||^2 + z * ||x||_1
        x=HandleProx(v,L);
        x_minus_s=x-s;
        r_sum = norm(x_minus_s)^2;
        f_new = HandleObjSmooth(x);
        l_sum = f_new - f_old - g_old'*x_minus_s;
%         if (r_sum <=1e-100)
%             flag=1; % this shows that, the gradient step makes little improvement
%             break;
%         end
        if(l_sum <= 0.5*r_sum*L)
            break;
        else
            L=2*L;
        end
    end
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
    beta=(alphap-1)/alpha; s = x + beta* xxp; xxp = x-xp;
    f_curr = f_new  + HandleObjNonSmooth(x);
    if(f_curr<f_best)
        f_best = f_curr;
        x_best = x;
    end
    his =[his;f_curr];
    if(flag==1),break;end
end




