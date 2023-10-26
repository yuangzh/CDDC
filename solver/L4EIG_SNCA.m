function [x,fobjs,ts] = L4EIG_SNCA(x,A,beta,stopCriteria,timeLimit,timeIntervel)
% min_x  x'Ax x'x + beta ||x||_4^4 - ||x||_2^2, s.t. -R <= x <= R
% x = Rescale_CD_SNCA(x,A,beta);

initt = clock;
last_rec_clock = initt;

n = size(A,1);
R = sqrt(beta/n);
% R = inf;
x = min(R,x);
x = max(-R,x);

% x = Rescale_CD_SNCA(x,A,beta);
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
    xAx1 = xAx-1;  a = A(i,i);
    
    %    c0 = xAx1*xtx  + beta*xxxx;
    c1 = xAx1*2*x(i) + 2*Ax(i)*xtx + 4*beta* xx(i)*x(i);
    c2 = xAx1  + 4*Ax(i)*x(i) + a*xtx  + 2*beta*xx(i)   + 4*beta* x(i)^2;
    c3 = 2*Ax(i) + 2*a*x(i) + 4*beta*x(i) ;
    c4 = a  + beta;
    %   c0 + c1*t + c2*t*t + c3*t*t*t + c4*t*t*t*t          -    HandleObj_t(t)
    
    lb = -R-x(i);  ub = R - x(i);
    [t_best] = ncvx_prox_poly4(c1,c2,c3,c4,lb,ub);
    x(i) = x(i)+ t_best;
    Ax = Ax + t_best*A(:,i);
    xx(i) = x(i)*x(i);
    xtx = xtx + x(i)*t_best*2 - t_best*t_best;
    xAx = xAx + 2*t_best*Ax(i)-t_best*t_best*A(i,i);
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
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


function x = Rescale_CD_SNCA(x,A,beta)
% rescale x, it is equvalent to: min_t F(t*x)
% xtx = x'*x;
% norm_x = sqrt(xtx);
% t = norm_x  / sqrt( 2 * ( x'*A*x*xtx + beta * norm(x,4)^4) );
% x =  t*x;

x =  x*norm(x)  / sqrt( 2 * ( x'*A*x*x'*x + beta * norm(x,4)^4) );
