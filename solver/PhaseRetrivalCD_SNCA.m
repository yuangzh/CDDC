function [x,fobjs,ts] = PhaseRetrivalCD_SNCA(Q,normG,x,stopCriteria,timeLimit,timeIntervel)
% min_{x} - 0.5 x'Ax, s.t. x \in [-1,+1]^n
% where A = normG*eye(size(G,1)) - Q;

initt = clock;
last_rec_clock = initt;

A = 1.01*normG*eye(size(Q,1)) - Q;
n = length(x);

theta = 1e-6;
% min_{x} 0.5 v'A v, s.t. x \in [-1,+1]^n
% min_eta 0.5*theta*eta^2 - 0.5 (x+eta*ei)'A (x+eta*ei), s.t. -1 <= x+eta*ei < = 1
% min_eta 0.5*(theta-A(i,i))*eta^2 - (Ax)_i eta , s.t. -1 <= x+eta*ei < = 1

HandleObj = @(x)PhaseRetrival_ComputeTrueObj(x,Q);

fobj = HandleObj(x); fobjs = [fobj]; ts = [etime(clock,initt)];

Ax = A*x;
fobj_old = fobj;
for iter = 1:1e10
    if(~mod(iter,2))
        i = mod(iter/2,n)+1;    % i = randperm(n,1);
    else
        grad =  ncvx_prox_quad(theta,-Ax,-1 - x, 1 - x);
        [~,i] = max(abs(grad));
    end
    
    t_best = ncvx_prox_quad(theta-A(i,i),-Ax(i),-1 - x(i),1 - x(i));
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
