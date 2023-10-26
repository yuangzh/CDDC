function [x,fobjs,ts] = OneNet_CD_SCA(x,G,y,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5||max(Gx,0) - max(y,0) ||_2^2
% min_x 0.5||max(Gx,0)||_2^2  + 0.5y'y - <max(Gx,0),|y|>
% min_x 0.5||max(Gx,0)||_2^2  + const - sum(max(A*x,0)), where A = diag(y)*G

initt = clock;
last_rec_clock = initt;

n = length(x);
y = max(0,y);
HandleObj = @(x)OneNet_ComputeTrueObj(x,G,y);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

A = bsxfun(@times,y,G); % A = diag(y)*G;
c = sum(G.*G,1)'+ 1e-6;
Gx = G*x;
Ax = A*x;
Gxplus = max(0,Gx);

fobj_old = fobj;
for iter =1:1e10
    
    i = mod(iter,n)+1;
    % min_alpha 0.5 c(i) alpha^2 + alpha (grad_f (i) - subgrad_g(i))
    % min_alpha 0.5 alpha^2 + alpha [grad_f(i)-subgrad_g(i)] / c(i)
    %     Index = find(Ax>0);subgrad_g_i = sum(A(Index,i));
    subgrad_g_i = sum( A(Ax>0,i) );
    grad_f_i = Gxplus'*G(:,i)  ;
    t_best = - (grad_f_i - subgrad_g_i) / c(i);
    
    x(i) = x(i) + t_best;
    Ax = Ax + t_best*A(:,i);
    Gx = Gx + t_best*G(:,i);
    Gxplus = max(Gx,0);
    
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = 0.5*norm(Gxplus-y)^2 ; % fobj = HandleObj(x);
        %         fprintf('CD-SCA iteration:%d, fobj:%.5f\n',iter,fobj);
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end

    if(~mod(iter,n))
        fobj_new = 0.5*norm(Gxplus-y)^2 ;
        rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
%         fprintf('CD-SCA iter:%d, fobj:%f, change:%f\n',iter,fobj_new,rel_change);
        if(rel_change < stopCriteria),break;end
        fobj_old = fobj_new;
    end
    
    
end

fobj = HandleObj(x);
% fprintf('CD-SCA iteration:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

% fprintf('\n\n\n\n\n\n\n\n\n');
