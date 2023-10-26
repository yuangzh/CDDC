function [x,fobjs,ts] = PhaseRetrivalGPM(Q,normQ,x,stopCriteria,timeLimit,timeIntervel)
% min_{x} 0.5 x'Qx, s.t. x \in {-1,+1}^n
% min_{x} -0.5 x'Ax, s.t. x \in {-1,+1}^n

initt = clock;
last_rec_clock = initt;

n = length(x);
A = normQ*eye(n) - Q;
HandleObj = @(x) PhaseRetrival_ComputeTrueObj(x,Q);

fobj = HandleObj(x);
fobjs = [fobj];
ts = [etime(clock,initt)];
fobj_old = fobj;
for iter = 1:1e10
    x = sign(A*x);
    
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

