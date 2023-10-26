function [x,fobjs,ts] = BinaryOpt_CD_SCA(G,y,x,lambda,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (sqrt(n) - ||x||_2), s.t. -1 <= x <= 1

initt = clock;
last_rec_clock = initt;

n = length(x);
x = max(-1,min(x,1));

HandleObj = @(x) 0.5*norm(G*x-y)^2 + lambda * ( sqrt(n) - norm(x) );

fobj = HandleObj(x);
fobjs = [fobj];
ts = [etime(clock,initt)];

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
    
    % min_t 0.5 ||G(x+tei)-y||_2^2  - lambda <sub_grad,tei>, s.t. -1 <= x+tei <= 1
    % min_t 0.5 ||c + t d ||_2^2  - lambda <sub_grad,tei>, s.t. -1 <= x+tei <= 1
    % min_t 0.5 d'd t^2 + <c,d> t  - lambda <sub_grad,tei>, s.t. -1 <= x+tei <= 1
    % min_t 0.5 t^2 + <c,d>/d'd t  - lambda <sub_grad,tei>, s.t. -1 <= x+tei <= 1
    % min_t 0.5 t^2 + <c,d>/d'd t - lambda / d'd <sub_grad,tei>, s.t. -1-x(i) <= t <= 1-x(i)
    t = (GGx(i)-Gy(i)-lambda*x(i)/sqrt(xx))/gg(i);
    t_best = max(-1-x(i),min(1-x(i),-t));
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