function [x,fobjs,ts] = OneNet_CD_SNCA(x,G,y,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5||max(Gx,0)-max(y,0)||_2^2
% min_x 0.5||max(Gx,0)||_2^2 + 0.5y'y - <max(Gx,0),|y|>

initt = clock;
last_rec_clock = initt;

y_plus = max(y,0);
n = length(x);
HandleObj = @(x)OneNet_ComputeTrueObj(x,G,y);
theta = 1e-6;

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

Gx = G*x;
fobj_old = 1e100;
for iter =1:10000000
    i = mod(iter,n)+1;
    % Note: If you successfully compile the function ncvx_prox_relu_whileloop.cpp, please use ncvx_prox_relu_mexC.m instead of ncvx_prox_relu.m
    t_best = ncvx_prox_relu(theta,0,G(:,i),Gx,y_plus);
    % t_best = ncvx_prox_relu_mexC(theta,0,G(:,i),Gx,y_plus);
    x(i) = x(i) + t_best;
    % reconstrct Gx in linear time
    Gx = Gx + t_best*G(:,i);
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = 0.5*norm(max(Gx,0)-y_plus)^2;
        %                 fprintf('CD-SNCA iteration:%d, fobj:%.5f\n',iter,fobj);
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end
    
    if(~mod(iter,n))
        fobj_new = 0.5*norm(max(Gx,0)-y_plus)^2;
        rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));
        %           fprintf('CD-SNCA iter:%d, fobj:%f, change:%f\n',iter,fobj_new,rel_change);
        if(rel_change < stopCriteria),break;end
        fobj_old = fobj_new;
    end
end

fobj = HandleObj(x);
% fprintf('CD-SNCA iteration:%d, fobj:%.5f, time: %f\n',iter,fobj,etime(clock,initt));
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

