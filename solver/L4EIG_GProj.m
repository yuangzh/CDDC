function [x,fobjs,ts] = L4EIG_GProj(x,A,beta,stopCriteria,timeLimit,timeIntervel)
% min_x x'Ax + beta ||x||_4^4, s.t. ||x|| = 1
% grad = 2*A*x + 4*beta*x.*x.*x;
% hess = 2*A + 12*beta*diag(x.*x)

initt = clock;
last_rec_clock = initt;

% normA = SpectralNorm(A);

HandleObj = @(x)L4EIG_ComputeTrueObj(x,A,beta);
HandleObjSmooth = @(x) ComputeObjSmooth(x,A,beta);
HandleProj = @(x) x/ norm(x);
x = HandleProj(x);
% fobjs = [];
% ts = [];
% fobj = HandleObj(x);
% fobjs = [fobjs;fobj];
% ts = [ts;etime(clock,initt)];

n=length(x); flag=0;xp=x;xxp=zeros(n,1);alpha=1; s=x;

fobj = L4EIG_ComputeTrueObj(x,A,beta);
fobjs = [fobj];
ts = [etime(clock,initt)];

fobj_old = fobj;
for iterStep=1:100000
    %     fprintf('*');
    [f_old,g_old] = HandleObjSmooth(s);
    x_old=x;
    max_in = 10000;L = 0.1;
    for in=1:max_in,
        v=s-g_old/L;
        % min_{x} 0.5 L ||x-(xt-g/L)||^2 + z * ||x||_1
        [x]=HandleProj(v);
        
        xs=x-s;
        r_sum= norm(xs)^2;
        f_new = HandleObjSmooth(x);
        l_sum = f_new - f_old - g_old'*xs;
        if (r_sum <=1e-100)
            flag=1; % this shows that, the gradient step makes little improvement
            break;
        end
        if(l_sum <= 0.5*r_sum*L)
            break;
        else
            %             L=max(2*L,l_sum/r_sum);
            L=2*L;
        end
        if(in==max_in)
            fprintf('warning! Lipschitz too large!');
        end
    end
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
    s=x + (alphap-1)/alpha* xxp;
    xxp=x-x_old;

    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
        fobjs  = [fobjs;fobj];
        ElasTime =  etime(cur_clock,initt);
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




function [fobj,grad] = ComputeObjSmooth(x,A,beta)
Ax = A*x;
xx = x.*x;
fobj = x'*Ax + beta*sum(xx.*xx);
grad = 2*Ax + beta*4*xx.*x;

