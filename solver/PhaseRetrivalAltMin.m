function [x,fobjs,ts] = PhaseRetrivalAltMin(G,Q,y,x,stopCriteria,timeLimit,timeIntervel)
% min_x min_s 0.5 ||Gs - xoy||_2^2


pinvG = pinv(full(G));

initt = clock;
last_rec_clock = initt;

HandleObj = @(x) PhaseRetrival_ComputeTrueObj(x,Q);
fobj = HandleObj(x);
fobjs = [fobj];
ts = [etime(clock,initt)];
fobj_old = fobj;

for iter = 1:200000
    
    % fix x, solve s
    s = pinvG*(x.*y);
    
    % fix s, solve x
    Gs = G*s;
    x = ProjBinary(Gs./y);
    
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



function [v] = ComputeOptimalxPhaseRetrival(G,s,y)
% 0.5 ||Gs - xoy||_2^2, x \in {-1,+1}
v = (G*s)./y;
v = v./abs(v);