function [x,fobjs,ts] = SparseOpt_CD_SNCA(G,y,x,lambda,k,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (||x||_1 - ||x||_{top_k})
% min_x 0.5x'G'Gx - <Gx,y> + lambda (||x||_1 - ||x||_{top_k})
% 0.5 t^2 + g_i/L t + lambda/L (||x+tei||_1 - ||x+tei||_{top_k})

initt = clock;
last_rec_clock = initt;

n = length(x);
HandleObj = @(x) 0.5*norm(G*x-y)^2 + lambda * ( norm(x,1) - topksum(x,k));

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

theta = 1e-6;

GG = G'*G;
gg =  diag(GG) + theta;
Gy = G'*y;
GGx = GG*x;
fobj_old = fobj;

for iter = 1:1e10
    subgrad_minus_g = - lambda*top_k_subgrad(x,k);
    if(~mod(iter,2))
        i = mod(iter/2,n)+1;
    else
        grad = prox_l1( x - ( GGx - Gy + subgrad_minus_g)./ gg, lambda./gg) - x;
        [~,i] = max(abs(grad));
    end
    
    % min_t 0.5 gg(i) t^2 +  (GGx(i) - Gy(i))* t + lambda ||x+tei||_1 - lambda || x+tei ||_{top-k}
    t_best =  ncvx_prox_l1topk(gg(i),GGx(i) - Gy(i),lambda,x,i,k);
    x(i) = x(i)+t_best;
    GGx = GGx + t_best*GG(:,i);
    
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








