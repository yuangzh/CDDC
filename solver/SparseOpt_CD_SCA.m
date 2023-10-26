function [x,fobjs,ts] = SparseOpt_CD_SCA(G,y,x,lambda,k,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (||x||_1 - ||x||_{top_k})
% min_x 0.5 ||Gx-y||_2^2  + lambda (||x||_1 - ||x||_{top_k})
% 0.5 x'GGx - <x,G'y> + lambda (||x||_1 - ||x||_{top_k})

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
    
    % min_t 0.5 (x+tei)'Q(x+tei) + <x+tei,p> + lambda || x+tei||_1 + <x+tei, subgrad_minus_g>
    % min_t 0.5 Q(i,i) t^2 + ((Qx)_i + p_i +subgrad_minus_g_i ) t+ lambda || x+tei||_1
    cof_b = GGx(i) - Gy(i) + subgrad_minus_g(i);
    % min_t 0.5 Q(i,i) t^2 + cof_b t+ lambda | x(i) + t |_1
    % x(i) + t = theta
    % min_theta 0.5 Q(i,i) (theta-x(i))^2 + cof_b (theta-x(i)) + lambda | theta |_1
    % min_theta 0.5 (theta-x(i))^2 + cof_b/Q(i,i) theta + lambda/Q(i,i) | theta |_1
    % min_theta 0.5 theta^2 - x(i)*theta + cof_b/Q(i,i) theta + lambda/Q(i,i) | theta |_1
    
    t_best = prox_l1(x(i)-cof_b/gg(i),lambda/gg(i)) - x(i);
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
