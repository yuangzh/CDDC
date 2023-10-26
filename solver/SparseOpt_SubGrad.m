function [x,fobjs,ts] = SparseOpt_SubGrad(G,y,x,lambda,k,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (||x||_1 - ||x||_{top_k})

initt = clock;
last_rec_clock = initt;

n = length(x);
HandleObj   = @(x) SparseOpt_ComputeTrueObj(x,G,y,lambda,k);
proj = @(x)max(min(x,1e6),-1e6);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

fobj_old = fobj;
for iter = 1:1e10
    
%     subgrad = top_k_subgrad(x,k);
    
    grad = G'*(G*x-y) + lambda*ComputeSubgrad(x,k);
    step = 0.01 / (iter);
    x_old = x;
    x = x - step*grad;
    x = proj(x);
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
    if(fobj_new >fobj_old)
        x = x_old;
    end
    rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
    if(rel_change < stopCriteria),break;end
    fobj_old = fobj_new;
    
    
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

function [subgrad] = ComputeSubgrad(x,k)
% Compute the subgradient of the function ||x||_1 - ||x||_{top_k}
n = length(x);
[~,ind]=sort(abs(x),'ascend');
I = ind(1:(n-k));
subgrad = zeros(n,1);
subgrad(I) = sign(x(I));


