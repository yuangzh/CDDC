function [x,fobjs,ts] = PhaseRetrivalGradProj(Q,normQ,x,stopCriteria,timeLimit,timeIntervel)
% min_{x} 0.5 x'Qx, s.t. x \in {-1,+1}^n
% x = min_x 0.5 normQ ||x-xt||_2^2 + <x-xt, grad f(xt)>, s.t. x \in {-1,+1}^n

initt = clock;
last_rec_clock = initt;

L = normQ;
n = length(x);
HandleObj = @(x)PhaseRetrival_ComputeTrueObj(x,Q);
HandleObjGP = @(x)ComputeObj(x,Q);

fobj = HandleObjGP(x);
fobjs = [fobj];
ts = [etime(clock,initt)];
fobj_old = fobj;
for iter = 1:1e10
    
    x_old = x;
    [f_old,g_old] = HandleObjGP(x);
    
    L = 0.01;
    for line_search = 1:100000
        v = x_old-g_old/L;
        [x] = ProjBinary(v);
        f_new = HandleObjGP(x);
        delta_x = x - x_old;
        l_sum = f_new - f_old - g_old'*delta_x;
        if(l_sum <= 0.5*L*norm(delta_x)^2)
            break;
        else
            L = L*1.1;
        end
    end
    
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


function [fobj,grad] = ComputeObj(x,Q)
Qx = Q*x;
grad = Qx;
fobj = 0.5* x'*Qx;



