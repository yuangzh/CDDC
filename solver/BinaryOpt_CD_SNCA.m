function [x,fobjs,ts] = BinaryOpt_CD_SNCA(G,y,x,lambda,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (sqrt(n) - ||x||_2), s.t. -1 <= x <= 1
% min_x 0.5 x'G'Gx - <Gx,y> + 0.5y'y + lambda (sqrt(n) - ||x||_2), s.t. -1 <= x <= 1
initt = clock;
last_rec_clock = initt;

n = length(x);
x = max(-1,min(x,1));

HandleObj = @(x) 0.5*norm(G*x-y)^2 + lambda * ( sqrt(n) - norm(x) );

fobj = HandleObj(x); fobjs = [fobj]; ts = [etime(clock,initt)];

theta = 1e-6;
GG = G'*G;
gg =  diag(GG) + theta;
Gy = G'*y;
GGx = GG*x;
xx = x'*x;
fobj_old = fobj;
for iter = 1:1e10
    
    if(~mod(iter,2))
        i = mod(iter/2,n)+1;
    else
        grad = max(-1-x,min(1-x, (-GGx + Gy + lambda/sqrt(xx)*x)./gg));
        [~,i] = max(abs(grad));
    end
    
    t_best = ncvx_prox_linfl2(gg(i),GGx(i) - Gy(i),lambda,x,i,-1-x(i),1-x(i));
    x(i) = x(i) + t_best;
    GGx = GGx + t_best*GG(:,i);
    xx = xx + 2*t_best*x(i) - t_best*t_best;
    
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