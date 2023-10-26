function [x,fobjs,ts] = BinaryOpt_PDCA(G,y,x,lambda,stopCriteria,timeLimit,timeIntervel)% min_x 0.5 ||Gx-y||_2^2  + lambda (sqrt(n) - ||x||_2), s.t. -1 <= x <= 1% min_x min_y 0.5 ||Gx-y||_2^2  + lambda (sqrt(n) - <x,y>), s.t. -1 <= x <= 1, ||y||_2^2 \leq ninitt = clock;last_rec_clock = initt;n = length(x);% L = SpectralNorm(G);x = max(-1,min(x,1));HandleObj = @(x)BinaryOpt_ComputeTrueObj(x,G,y,lambda,n);fobj = HandleObj(x);fobjs = [fobj];ts = [etime(clock,initt)];s = x; xxp = zeros(n,1); alpha = 1;fobj_old = fobj;for iter = 1:1e10            HandleObjSmooth = @(x) ComputeSubProbObjGrad(x,G,y,lambda,n);    [f_old,g_old] = HandleObjSmooth(s);        x_old = x;    L = 0.1;    for i_linesearch = 1:1000,        x = s - g_old/L;        x = max(-1,min(x,1));        xs = x-s;        r_sum = norm(xs)^2;        f_new = HandleObjSmooth(x);        l_sum = f_new - f_old - g_old'*xs;        if(l_sum <= 0.5*r_sum*L),break;        else L=2*L; end    end    alphap=alpha; alpha=(1+sqrt(4*alpha*alpha +1))/2;    s = x + (alphap-1)/alpha*xxp; xxp=x-x_old;        % min_x 0.5 ||Gx-y||_2^2  + <x,subgrad>, s.t. -1 <= x <= 1    %     subgrad = lambda*x / norm(x);    %     grad = G'*(G*x-y) - subgrad;    %     x = min(1,max(x - grad/L,-1));            cur_clock = clock;    if(etime(cur_clock,last_rec_clock) > timeIntervel)        fobj = HandleObj(x);        ElasTime =  etime(cur_clock,initt);        fobjs  = [fobjs;fobj];        ts = [ts;ElasTime];        last_rec_clock = cur_clock;        if ElasTime > timeLimit            break;        end    end        fobj_new = HandleObj(x);    rel_change = abs(fobj_old - fobj_new) / (1+abs(fobj_old));    if(rel_change < stopCriteria),break;end    fobj_old = fobj_new;endfobj = HandleObj(x);% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);fobjs  = [fobjs;fobj];ts = [ts;etime(clock,initt)];function [fobj,grad_dc] = ComputeSubProbObjGrad(x,G,y,lambda,n)diff = G*x-y;fobj = 0.5*norm(diff)^2 + lambda*(sqrt(n) - norm(x));grad_dc = G'*diff - lambda*x/norm(x);