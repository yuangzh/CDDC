function [x,fobjs,ts] = L1PCA_MSCR(G,x,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 x'*x - ||Gx||_1
% min_x 0.5 x'*x - <x,subgrad>

initt = clock;
last_rec_clock = initt;

HandleObj = @(x)L1PCA_ComputeTrueObj(x,G);

fobjs = []; ts = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts  = [ts;etime(clock,initt)];
fobj_old = fobj;
for iter =1:1e10
    
    subgrad = G'*sign(G*x);
    x = subgrad;
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
%         fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
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
