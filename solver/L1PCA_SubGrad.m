function [x,fobjs,ts] = L1PCA_SubGrad(G,x,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 x'*x - ||Gx||_1
% min_x 0.5 x'*x - <x,subgrad>

initClock = clock;
lastClock = initClock;

HandleObj = @(x)0.5*x'*x - norm(G*x,1);

fobj = HandleObj(x);
fobjs = [fobj];
ts = [etime(clock,initClock)];
fobj_old = fobj;
for iter =1:10000
    
    % update x
    subGrad = x - G'*sign(G*x);
    stepSize = 1 / (iter);
    x = x - stepSize*subGrad;
    
    currClock = clock;
    if(etime(currClock,lastClock) > timeIntervel)
        fobj = HandleObj(x);
        timeElapsed  =  etime(currClock,initClock);
        fobjs  = [fobjs;fobj];
        ts = [ts;timeElapsed ];
        lastClock = currClock;
        if timeElapsed  > timeLimit
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
ts = [ts;etime(clock,initClock)];