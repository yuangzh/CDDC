function [x,fobjs,ts] = L4EIG_PDCA(x,A,beta,stopCriteria,timeLimit,timeIntervel)
% min_x  (x'Ax x'x + beta ||x||_4^4) -  (||x||_2^2), s.t. -R <= x <= R

initt = clock;
last_rec_clock = initt;

n = size(A,1);
min_eig_A = 0;
A = A + min_eig_A*eye(n);


% x = Rescale_CD_SNCA(x,A,beta);
r1 = 1 / (min_eig_A + beta/n);
r1_sqrt = sqrt(r1);

L = 1;
HandleObj = @(x)L4EIG_ComputeTrueObj(x,A,beta);
HandleObjSmooth = @(x) ComputeObjSmooth2(x,A,beta);
HandleProj = @(x) min(r1_sqrt,max(-r1_sqrt,x));
x = HandleProj(x);

flag=0; xp=x; xxp=zeros(n,1); alpha=1; s=x;

fobj = L4EIG_ComputeTrueObj(x,A,beta);
fobjs = [fobj];
ts = [etime(clock,initt)];

fobj_old = fobj;
for iterStep=1:100000
    %     fprintf('*');
    [f_old,g_old] = HandleObjSmooth(s);
    x_old=x;
    max_in = 10000;
    for in=1:max_in,
        v=s-g_old/L;
        % min_{x} 0.5 L ||x-(xt-g/L)||^2 + z * ||x||_1
        [x]=HandleProj(v);
        xs=x-s;
        r_sum=norm(xs)^2;
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
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2; alpha = 1;
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


function [fobj,grad] = ComputeObjSmooth2(x,A,beta)
Ax = A*x;
xx = x.*x;
xtx = sum(xx);
xAx = x'*Ax;
fobj = xAx*xtx + beta*norm(x,4)^4 - x'*x;
grad = 2*Ax*xtx + 2*xAx*x + 4*beta*x.*xx - 2*x;


function x = Rescale_CD_SNCA(x,A,beta)
% rescale x, it is equvalent to: min_t F(t*x)
% xtx = x'*x;
% norm_x = sqrt(xtx);
% t = norm_x  / sqrt( 2 * ( x'*A*x*xtx + beta * norm(x,4)^4) );
% x =  t*x;

x =  x*norm(x)  / sqrt( 2 * ( x'*A*x*x'*x + beta * norm(x,4)^4) );

