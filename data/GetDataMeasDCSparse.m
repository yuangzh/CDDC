function [y] = GetDataDCSparse_y(G)
rand('seed',0);
randn('seed',0);
[m,n] = size(G);
totnnz = 100;
x = generate_x(n,totnnz);
b = G*x;
y = b + 0.1*norm(b)*randn(m,1);
 

function x = generate_x(n,k)
p = randperm(n);
x = zeros(n,1);
x(p(1:k)) = randn(k,1);