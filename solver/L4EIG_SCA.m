function [x,fobjs,ts] = L4EIG_SCA(x,A,beta,stopCriteria,timeLimit,timeIntervel)
% min_x  sigma x'x - (sigma x'x - x'Ax x'x - beta ||x||_4^4 + ||x||_2^2), s.t. -R <= x <= R

initt = clock;
last_rec_clock = initt;

n = size(A,1);
min_eig_A = 0;
A = A + min_eig_A*eye(n);


% x = Rescale_CD_SNCA(x,A,beta);
r1 = 1 / (min_eig_A + beta/n);
r1_sqrt = sqrt(r1);
r2 = 1 / ( min_eig_A + beta / (n*n));
sigma = 6*beta*r1 + 6*sqrt(SpectralNorm(A))*r2  - 1;

x = min(r1_sqrt,x);
x = max(-r1_sqrt,x);

HandleObj = @(x)L4EIG_ComputeTrueObj(x,A,beta);
HandleObj2 = @(x)x'*A*x*x'*x + beta*sum(x.*x.*x.*x) - x'*x;

fobj = L4EIG_ComputeTrueObj(x,A,beta);
fobjs = [fobj];
ts = [etime(clock,initt)];

Ax = A*x;
xx = x.*x;
xtx = sum(xx);
xAx = x'*Ax;
fobj_old = fobj;
for iter = 1:3000000000

    if(~mod(iter,2))
        i = mod(iter/2,n)+1;
    else
        grad = 2*Ax*xtx + 2*xAx*x + 4*beta*x.*xx - 2*x;
        [~,i] = max(abs(grad));
    end

    grad_g_i = 2*sigma*x(i) - 2*Ax(i)*xtx - 2*xAx*x(i) - 4*beta*x(i).*xx(i) + 2*x(i);
    % min_t theta ||x + t *ei ||_2^2 - <t*ei,grad>
    % min_t theta t*t + theta*2*x(i)*t) - grad(i)*t
    lb = -r1_sqrt-x(i);
    ub = r1_sqrt - x(i);
    t_best = min_quad(sigma,sigma*2*x(i)-grad_g_i,lb,ub);
    x(i) = x(i) + t_best;
    
    Ax = Ax + t_best*A(:,i);
    xx(i) = x(i)*x(i);
    xtx = xtx + x(i)*t_best*2 - t_best*t_best;
    xAx = xAx + 2*t_best*Ax(i)-t_best*t_best*A(i,i);
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
        fobjs  = [fobjs;fobj];
        ElasTime =  etime(cur_clock,initt);
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end
    
    if(~mod(iter,2*n))
        fobj_new = HandleObj(x);
        rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
        if(rel_change < stopCriteria),break;end
        fobj_old = fobj_new;
    end
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

function t = min_quad(a,b,lb,ub)
% min_t a*t*t + b*t, s.t. lb<=t<= ub
fx = @(t)a*t*t+b*t;
t1 = -b / (2*a);
t2 = lb;
t3 = ub;
f1 = fx(t1);
f2 = fx(t2);
f3 = fx(t3);
if(f1<=f2 && f1<= f3)
    t = t1;
elseif(f2<=f1 && f2<=f3)
    t = t2;
else
    t = t3;
end


function x = Rescale_CD_SNCA(x,A,beta)
% rescale x, it is equvalent to: min_t F(t*x)
% xtx = x'*x;
% norm_x = sqrt(xtx);
% t = norm_x  / sqrt( 2 * ( x'*A*x*xtx + beta * norm(x,4)^4) );
% x =  t*x;

x =  x*norm(x)  / sqrt( 2 * ( x'*A*x*x'*x + beta * norm(x,4)^4) );

