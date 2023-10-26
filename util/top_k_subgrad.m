function subgrad = top_k_subgrad(x,k)
% Compute the subgradient of the top k function (in absolute value)
% min_v  -<x,v>, s.t. v \in {-1,+1}, ||v'||_1 \leq k
% subgrad = v = sign(x)w;
% min_w <abs(x),w>, s.t. w \in {0,1}, sum(w)=k
n = length(x);
[~,ind]=sort(abs(x),'descend');
w = zeros(n,1);
w(ind(1:k))=1;
subgrad = sign1(x).*w;

function [v] = sign1(v)
% min_{v} || v - a||_2^2, s.t. v \in{-1,+1}^n
absv = abs(v);
v = v./absv;
v(absv==0) = 1;

