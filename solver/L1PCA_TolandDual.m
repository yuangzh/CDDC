function [x,fobjs,ts] = L1PCA_TolandDual(G,x,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 x'*x - ||Gx||_1
% min_x min_y 0.5 x'*x - <y,Gx>, s.t. -1<=y<=1
% min_y min_x 0.5 x'*x - <y,Gx>, s.t. -1<=y<=1
% grad = x - G'y = 0  => x = G'y
% min_y -0.5 y'*G*G'*y, s.t. -1<=y<=1
initt = clock;
last_rec_clock = initt;

y = sign(G*x);
HandleObj = @(x)0.5*x'*x - norm(G*x,1);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];
fobj_old = fobj;
for iter = 1:1e10
    
    y = sign(G*G'*y);
    x = G'*y;
    
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

x = G'*y;
fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];