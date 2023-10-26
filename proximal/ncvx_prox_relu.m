%{

% Note: You can use the following code to validate the correctness of the algorithm

function test_ncvx_prox_relu
clc;clear all;close all;
mex ncvx_prox_relu_whileloop.cpp

for iter = 1:10000000
    
    rand('seed',iter);
    randn('seed',iter);
    n = 15;
    alpha = rand(1)*100*rand(1);
    beta = randn(1)*100*rand(1);
    a = randn(n,1)*100*rand(1).*max(0,randn(n,1));
    b = randn(n,1)*100*rand(1).*max(0,randn(n,1));
    w = rand(n,1)*100*rand(1).*max(0,randn(n,1))*0;
    
    fun = @(t)0.5*alpha*t^2 + beta*t + 0.5*norm(max(0,t*a + b))^2  - w'*max(0,t*a + b);
    t1 = fminbnd(fun,-100000,100000);
    t2 = ncvx_prox_relu___(alpha,beta,a,b,w);
    f1 = fun(t1);
    f2 = fun(t2);
    fprintf('iter:%.6d, fobj:%.1e %.1e, diff: %.5e\n',iter,f1,f2,f1-f2);
    if(f2>f1+1e-5*abs(mean([f1;f2])))
        x1
        x2
        dddd
    end
end

%}


function best_t = ncvx_prox_relu___(alpha,beta,a,b,w)
% min_t 0.5*alpha*t^2+beta*t+0.5*||max(0,t*a+b)||_2^2-w'*max(0,t*a+b)

I=find(a==0); a(I)=[]; b(I)=[]; w(I)=[];  a = full(a);
I1=find(a>0); a1=a(I1); b1=b(I1); w1=w(I1);
I2=find(a<0); a2=a(I2); b2=b(I2); w2=w(I2);

[~,o1]=sort(b1./a1,'ascend'); a1=a1(o1); b1=b1(o1); w1=w1(o1);
[~,o2]=sort(b2./a2,'ascend'); a2=a2(o2); b2=b2(o2); w2=w2(o2);

q1=a1.*a1; p1=a1.*(w1-b1); q2=a2.*a2; p2=a2.*(w2-b2);
i1=length(a1); i2=length(a2); z=[-1e20;sort(-b./a,'ascend');1e20];

if(i2==0), r2=0; s2=0; else r2=sum(q2); s2=sum(p2);  end
r1=0; s1=0; fobjs=[]; ts=[];

len_a = length(a);
best_t = ncvx_prox_relu_whileloop(alpha,beta,a,b,w,z,q1,q2,p1,p2,a1,a2,b1,b2,len_a,i1,i2,r1,r2,s1,s2);

return;

hf=@(t)0.5*alpha*t^2+beta*t+0.5*norm(max(0,t*a+b))^2-w'*max(0,t*a+b);
for i = 1:(length(a)+1)
    tt=(z(i)+z(i+1))/2;
    while ((i1>0)&&(tt*a1(i1)+b1(i1)>0))
        r1=r1+q1(i1); s1=s1+p1(i1); i1=i1-1;  end
    while ((i2>0)&&(tt*a2(i2)+b2(i2)<=0))
        r2=r2-q2(i2); s2=s2-p2(i2); i2=i2-1;  end
    tt=(s1+s2-beta)/(r1+r2+alpha);
    if((tt>=z(i)-1e-16)&&(tt<=z(i+1)+1e-16))
        fobjs=[fobjs;hf(tt)]; ts=[ts;tt];  end
end
[~,j] = min(fobjs);best_t = ts(j);





