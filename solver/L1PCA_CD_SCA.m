function [x,fobjs,ts] = L1PCA_CD_SCA(G,x,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 x'x - ||Gx||_1

initt = clock;
last_rec_clock = initt;

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
    subgrad_i = sign(Gx)'*G(:,i);
    
    %     ttt = min_t 0.5 (1+theta) t^2 + <x-grad,tei>
    %     ttt = min_t 0.5 t^2 +(x_i-grad_i)*t / (1+theta)
    t_best = (subgrad_i - x(i))/(1+theta);
    
    x(i) = x(i) + t_best;
    xx = xx + 2*t_best*x(i) - t_best*t_best;
    Gx = Gx + t_best*G(:,i);
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = 0.5*xx - norm(Gx,1);
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
        rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
        if(rel_change < stopCriteria),break;end
        fobj_old = fobj_new;
    end
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];
