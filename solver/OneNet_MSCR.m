function [x,fobjs,ts] = OneNet_MSCR(x,G,y,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5||max(Gx,0) - max(y,0) ||_2^2
% min_x 0.5||max(Gx,0)||_2^2  + 0.5y'y - <max(Gx,0),|y|>
% min_x 0.5||max(Gx,0)||_2^2  + const - sum(max(A*x,0)), where A = diag(y)*G

initt = clock;
last_rec_clock = initt;

const = 0.5*norm(y)^2;
y = max(y,0);
A = bsxfun(@times,y,G); % A = diag(y)*G;
% L = SpectralNorm(G);
HandleObj = @(x)OneNet_ComputeTrueObj(x,G,y);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

for iter = 1:1e10
    
    % Compute the subgradient of g(x)=sum(max(A*x,0))
    subgrad_g = sum(A((A*x)>0,:),1)';
    
    % Solve the following convex subproblem using APG
    % 0.5*|| max(Gx,0) ||_2^2 + const - <x - x^t,subgrad_g>
    HandleObjSmooth = @(x) ComputeSubProbObjGrad(x,G,const,subgrad_g);
    HandleObjNonSmooth = @(x)ComputeSubProbNonSmoothObj(x);
    HandleProx = @(a,theta)ComputeSubProbProxMap(a,theta);
    maxiter = 10;
    x_old = x;
    x = APG(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,maxiter);

    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
%         fprintf('MSCR iteration:%d, fobj:%.5f, diff:%f \n',iter,fobj,norm(x-x_old,1));
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end
    f_old = HandleObj(x_old);
    f_new = HandleObj(x);
    rel_change = abs(f_old - f_new) / (1+abs(f_old));
    if(rel_change < stopCriteria),break;end

end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];
        

function [fobj] = ComputeSubProbNonSmoothObj(x)
% For this example, there is no nonsmooth term
fobj = 0;

function [x] = ComputeSubProbProxMap(a,theta)
% For this example, g(x) = 0
% 0.5 theta ||x - a||^2 + g(x)
x = a;

function [fobj,grad] = ComputeSubProbObjGrad(x,G,const,subgrad_g)
Gx = G*x;
zo = Gx>0;
fobj = 0.5*norm(max(Gx,0))^2  + const - x'*subgrad_g;
grad = G'*(Gx.*zo) - subgrad_g;
