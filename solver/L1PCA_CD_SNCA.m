function [x,fobjs,ts] = L1PCA_CD_SNCA(G,x,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 x'x - ||Gx||_1

initt = clock;
last_rec_clock = initt;

% Unbaised = @(x)(norm(G*x,1)/(x'*x))*x;
% x = Unbaised(x);
        
n = length(x);
HandleObj = @(x)0.5*x'*x - norm(G*x,1);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

Gx = G*x;
xx = x'*x;
theta = 1e-6;
fobj_old = fobj;
for iter = 1:1000000
    i = mod(iter,n)+1;
    cof_alpha = 1+theta;
    cof_beta = x(i);
    cof_lambda = 1;
    cof_g = G(:,i);
    cof_d = Gx;
    t_best = ncvx_prox_l1(cof_alpha,cof_beta,cof_lambda,cof_g,cof_d);
    x(i) = x(i) + t_best;
    xx = xx + 2*t_best*x(i) - t_best*t_best;
    Gx = Gx + t_best*G(:,i);
 
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj =  0.5*xx - norm(Gx,1);
%         fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end
 
    if(~mod(iter,n))
        fobj_new = 0.5*xx - norm(Gx,1);
%          fprintf('iter:%d, fobj:%.5f\n',iter,fobj_new);
        rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
        if(rel_change < stopCriteria),break;end
        fobj_old = fobj_new;
    end
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

        
