function [x,fobjs,ts] = SparseOpt_PDCA(G,y,x,lambda,k,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5 ||Gx-y||_2^2  + lambda (||x||_1 - ||x||_{top_k})
% min_x min_y 0.5 ||Gx-y||_2^2  + lambda (||x||_1 - <x,y>), s.t. -1<=y<=1, ||y||_1 \leq k

initt = clock;
last_rec_clock = initt;

% L = SpectralNorm(G);
HandleObj   = @(x) SparseOpt_ComputeTrueObj(x,G,y,lambda,k);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

n = length(x);
s = x; xxp = zeros(n,1); alpha = 1;

fobj_old = fobj;
for iter = 1:1e10
    
    HandleObjSmooth = @(x) ComputeSubProbObjGrad(x,G,y,lambda,k);
    [f_old,g_old] = HandleObjSmooth(s);
    
    x_old = x;
    L = 1;
    for i_linesearch = 1:1000,
        x = prox_l1(s - g_old/L,lambda/L); xs = x-s;
        r_sum = norm(xs)^2;
        f_new = HandleObjSmooth(x);
        l_sum = f_new - f_old - g_old'*xs;
        if(l_sum <= 0.5*r_sum*L),break;
        else L=2*L; end
    end
    alphap=alpha; alpha=(1+sqrt(4*alpha*alpha +1))/2; alpha =1; % we set to alpha=1 since it is more stable
    s = x + (alphap-1)/alpha*xxp; xxp=x-x_old;
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
        %         fprintf('PDCA iteration:%d, fobj:%.5f\n',iter,fobj);
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


function [fobj,grad_dc] = ComputeSubProbObjGrad(x,G,y,lambda,k)
diff = G*x-y;
[ftopksum,subgrad] = topksum(x,k);
fobj = 0.5*norm(diff)^2 + lambda*(norm(x,1) - ftopksum);
grad_dc = G'*diff - lambda*subgrad;


