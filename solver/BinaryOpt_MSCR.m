function [x,fobjs,ts] = BinaryOpt_MSCR(G,y,x,lambda,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (sqrt(n) - ||x||_2), s.t. -1 <= x <= 1
% min_x min_y 0.5 ||Gx-y||_2^2  + lambda (sqrt(n) - <x,y>), s.t. -1 <= x <= 1, ||y||_2^2 \leq n
initt = clock;
last_rec_clock = initt;

n = length(x);
% L = SpectralNorm(G);
x = max(-1,min(x,1));

HandleObj = @(x)BinaryOpt_ComputeTrueObj(x,G,y,lambda,n);

fobj = HandleObj(x);
fobjs = [fobj];
ts = [etime(clock,initt)];
fobj_old = fobj;
for iter = 1:1e10
    
    subgrad_g = lambda*x / norm(x);
    
    HandleObjSmooth = @(x) ComputeSubProbObjGrad(x,G,y,subgrad_g);
    HandleObjNonSmooth = @(x)0;
    HandleProx = @(a,theta)ComputeSubProbProxMap(a,theta);
    maxiter = 10;
    x = APG(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,maxiter);
    
    % min_x 0.5 ||Gx-y||_2^2  + <x,subgrad_g>, s.t. -1 <= x <= 1
    % for in = 1:10
    %   grad = G'*(G*x-y) - subgrad_g;
    %   x_old = x;
    %   x = min(1,max(x - grad/L,-1));
    % end
    
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




function [fobj,grad] = ComputeSubProbObjGrad(x,G,y,subgrad_g)
diff = G*x-y;
fobj = 0.5*norm(diff)^2 - x'*subgrad_g;
grad = G'*diff - subgrad_g;

function [x] = ComputeSubProbProxMap(a,theta)
% 0.5 theta ||x - a||^2, s.t. -1 <= x <= 1
x = min(1,max(a,-1));