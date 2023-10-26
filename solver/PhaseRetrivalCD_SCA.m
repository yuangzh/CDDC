function [x,fobjs,ts] = PhaseRetrivalCD_SCA(Q,normG,x,stopCriteria,timeLimit,timeIntervel)
% min_{x} - 0.5 x'Ax, s.t. x \in [-1,+1]^n
% where A = normG*I - Q;

initt = clock;
last_rec_clock = initt;

n = length(x);
A = normG*eye(n) - Q;
HandleObj = @(x)PhaseRetrival_ComputeTrueObj(x,Q);

fobj = HandleObj(x); fobjs = [fobj]; ts = [etime(clock,initt)];

theta = 1e-5;
Ax = A*x;
fobj_old = fobj;

for iter = 1:1e10
    if(~mod(iter,2))
        i = mod(iter/2,n)+1;
    else
        grad =  ncvx_prox_quad(theta,-Ax,-1 - x, 1 - x);
        [~,i] = max(abs(grad));
    end
    
    t_best = ncvx_prox_quad(theta,-Ax(i),-1 - x(i), 1 - x(i));
    x(i) = x(i) + t_best;
    Ax = Ax + t_best*A(:,i);
    
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