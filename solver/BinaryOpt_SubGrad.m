function [x,fobjs,ts] = BinaryOpt_SubGrad(A,b,x,lambda,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Ax-b||_2^2  + lambda (sqrt(n) - ||x||_2), s.t. -1 <= x <= 1
% min_x min_y 0.5 ||Ax-b||_2^2  + lambda (sqrt(n) - <x,y>), s.t. -1 <= x <= 1, ||y||_2^2 \leq n
initt = clock;
last_rec_clock = initt;

n = length(x);
x = max(-1,min(x,1));

HandleObj = @(x)BinaryOpt_ComputeTrueObj(x,A,b,lambda,n);


fobj = HandleObj(x);
fobjs = [fobj];
ts = [etime(clock,initt)];
fobj_old = fobj;
for iter = 1:1e10
    
    grad = A'*(A*x-b) - lambda*x/norm(x);
    step = 1 / (iter);
    x = x - step*grad;
    x = max(-1,min(x,1));
    
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