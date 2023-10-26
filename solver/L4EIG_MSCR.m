function [x,fobjs,ts] = L4EIG_MSCR(x,A,beta,stopCriteria,timeLimit,timeIntervel)
% min_x  (x'Ax x'x + beta ||x||_4^4) -  (||x||_2^2), s.t. -r1_sqrt <= x <= r1_sqrt

initt = clock;
last_rec_clock = initt;

n = length(x);
min_eig_A = 0;
A = A + min_eig_A*eye(n);

HandleObj = @(x)L4EIG_ComputeTrueObj(x,A,beta);
% x = Rescale_CD_SNCA(x,A,beta);
r1 = 1 / (min_eig_A + beta/n);
r1_sqrt = sqrt(r1);

HandleProj = @(x) min(r1_sqrt,max(-r1_sqrt,x));
x = HandleProj(x);

fobj = L4EIG_ComputeTrueObj(x,A,beta);
fobjs = [fobj];
ts = [etime(clock,initt)];
fobj_old = fobj;
for iterStep=1:100000
    
    grad_g = 2*x;
    % Solve the following convex subproblem using APG
    % 0.5*|| max(Gx,0) ||_2^2 + const - <x-x^t,subgrad_g>
    HandleObjSmooth =  @(x) ComputeObjSmooth3(x,A,beta,grad_g);
    HandleObjNonSmooth = @(x)0;
    HandleProx = @(a,theta)ComputeSubProbProxMap(a,theta,r1_sqrt);
    maxiter = 10;
    x = APG(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,maxiter);
    
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
    
    fobj_new = HandleObj(x);
    rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
    if(rel_change < stopCriteria),break;end
    fobj_old = fobj_new;
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

function [x] = ComputeSubProbProxMap(a,theta,r1_sqrt)
% 0.5 theta ||x - a||^2, s.t. -r1_sqrt<=x<=r1_sqrt
x = max(min(a,r1_sqrt),-r1_sqrt);



function [fobj,grad] = ComputeObjSmooth2(x,A,beta)
Ax = A*x;
xx = x.*x;
xtx = sum(xx);
xAx = x'*Ax;
fobj = xAx*xtx + beta*norm(x,4)^4 - x'*x;
grad = 2*Ax*xtx + 2*xAx*x + 4*beta*x.*xx - 2*x;



function [fobj,grad] = ComputeObjSmooth3(x,A,beta,grad_g)
Ax = A*x;
xx = x.*x;
xtx = sum(xx);
xAx = x'*Ax;
fobj = xAx*xtx + beta*norm(x,4)^4 - x'*grad_g;
grad = 2*Ax*xtx + 2*xAx*x + 4*beta*x.*xx - grad_g;


function x = Rescale_CD_SNCA(x,A,beta)
% rescale x, it is equvalent to: min_t F(t*x)
% xtx = x'*x;
% norm_x = sqrt(xtx);
% t = norm_x  / sqrt( 2 * ( x'*A*x*xtx + beta * norm(x,4)^4) );
% x =  t*x;

x =  x*norm(x)  / sqrt( 2 * ( x'*A*x*x'*x + beta * norm(x,4)^4) );

