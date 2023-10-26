function [x,fobjs,ts] = SparseOpt_MSCR(G,y,x,lambda,k,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-b||_2^2  + lambda (||x||_1 - ||x||_{top_k})
% min_x min_y 0.5 ||Gx-b||_2^2  + lambda (||x||_1 - <x,y>), s.t. -1<=y<=1, ||y||_1 \leq k

initt = clock;
last_rec_clock = initt;

% L = SpectralNorm(G);
HandleObj = @(x)SparseOpt_ComputeTrueObj(x,G,y,lambda,k);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

fobj_old = fobj;

for iter = 1:1e10
    subgrad_g = lambda*top_k_subgrad(x,k);
    
    % min_x 0.5 ||Gx-y||_2^2 - x'subgrad_g  + lambda ||x||_1
    HandleObjSmooth = @(x) ComputeSubProbObjGrad(x,G,y,subgrad_g);
    HandleObjNonSmooth = @(x)lambda*norm(x,1);
    HandleProx = @(a,theta)ComputeSubProbProxMap(a,theta,lambda);
    maxiter = 10;
    x = APG(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,maxiter);
    
    %     for it = 1:10
    %         grad = G'*(G*x-y) - subgrad_g;
    %         x = prox_l1(x-grad/L,lambda/L);
    %     end
    
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

function [x] = ComputeSubProbProxMap(a,theta,lambda)
% 0.5 theta ||x - a||^2 + lambda ||x||_1
[x] = prox_l1(a,lambda/theta);



