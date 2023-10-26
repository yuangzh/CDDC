function [f,subgrad] = topksum(x,k)
n = length(x);
absx = abs(x);
[sort_absx,ind]=sort(absx,'descend');
% f = sum(absx(ind(1:k)));
f = sum(sort_absx(1:k));
w = zeros(n,1);
w(ind(1:k))=1;
subgrad = sign1(x).*w;



function [v] = sign1(v)
% min_{v} || v - a||_2^2, s.t. v \in{-1,+1}^n
absv = abs(v);
v = v./absv;
v(absv==0) = 1;


