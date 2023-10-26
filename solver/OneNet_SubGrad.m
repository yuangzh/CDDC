function [x,fobjs,ts] = OneNet_SubGrad(x,G,y,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5*||max(Gx,0) - |y| ||_2^2

initt = clock;
last_rec_clock = initt;

y = max(y,0);
HandleObj = @(x)OneNet_ComputeObjSubGrad(x,G,y);
proj = @(x)max(min(x,1e6),-1e6);
ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

for iter = 1:10000000
    [fobj,grad] = HandleObj(x);
%     if(iter==1)
%         scale = 1 / norm(grad(:));
%     end
    step = 1 / iter;
    x_old = x;
    x = x - step*grad;
    x = proj(x);
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
%         fprintf('SubGrad iteration:%d, fobj:%.5f\n',iter,fobj);
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end
    f_old = HandleObj(x_old);
    f_new = HandleObj(x);
    rel_change = abs(f_old - f_new) / (1+abs(f_old));
    if(rel_change < stopCriteria),break;end
%     if(norm(step*grad) <accuracy),break;end
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];


function [fobj,grad] = OneNet_ComputeObjSubGrad(x,G,y)
Gx = G*x;
Gx1 = max(0,Gx);
fobj = 0.5*norm(Gx1 - y)^2 ;
grad =  G'*(Gx1 - y) ;
