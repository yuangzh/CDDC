function [x,fobjs,ts] = OneNet_PDCA(x,G,y,stopCriteria,timeLimit,timeIntervel)
% min_x 0.5||max(Gx,0) - max(y,0) ||_2^2
% min_x 0.5||max(Gx,0)||_2^2  + 0.5y'y - <max(Gx,0),|y|>
% min_x 0.5||max(Gx,0)||_2^2  + const - sum(max(A*x,0)), where A = diag(y)*G

initt = clock;
last_rec_clock = initt;

n = length(x);
const = 0.5*norm(y)^2;
y = max(y,0);
A = bsxfun(@times,y,G); % A = diag(y)*G;
L = SpectralNorm(G);

HandleObj = @(x)OneNet_ComputeTrueObj(x,G,y);

ts = []; fobjs = [];
fobj = HandleObj(x);
fobjs = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

s = x; flag=0;

for iter = 1:1e10
    
    % Compute the subgradient of g(x)=sum(max(A*x,0))
    HandleObjSmooth = @(x) ComputeSubProbObjGrad(x,G,A,const);
    [f_old,g_old] = HandleObjSmooth(s);
    
    L = 1;
    for i_linesearch = 1:1000000
        x = s - g_old/L;
        xs = x-s;
        r_sum = norm(xs)^2;
        if (r_sum <=1e-10)
            flag=1;break;
        end
        f_new = HandleObjSmooth(x);
        l_sum = f_new - f_old - g_old'*xs;
        if(l_sum <= 0.5*r_sum*L),break;
        else L=2*L; end
    end
%     if(flag==1),break;end
    s = x ; 
    
    cur_clock = clock;
    if(etime(cur_clock,last_rec_clock) > timeIntervel)
        fobj = HandleObj(x);
        %         fprintf('PDCA iteration:%d, fobj:%.5f, diff: %f\n',iter,fobj,norm(xxp,1));
        ElasTime =  etime(cur_clock,initt);
        fobjs  = [fobjs;fobj];
        ts = [ts;ElasTime];
        last_rec_clock = cur_clock;
        if ElasTime > timeLimit
            break;
        end
    end
    rel_change = abs(f_old - f_new) / (1+abs(f_old));
    if(rel_change < stopCriteria),break;end
    %     if( norm(xxp,1) < accuracy),break;end
end

fobj = HandleObj(x);
% fprintf('iter:%d, fobj:%.5f\n',iter,fobj);
fobjs  = [fobjs;fobj];
ts = [ts;etime(clock,initt)];

function [fobj,grad_dc] = ComputeSubProbObjGrad(x,G,A,const)
Gx = G*x;
zo = Gx>0;
Ax = A*x;
fobj = 0.5*norm(max(Gx,0))^2  - sum(max(Ax,0)) + const;
subgrad_g = sum(A((Ax)>0,:),1)'; % subgradient of the convex function
grad_dc = G'*(Gx.*zo) - subgrad_g;
