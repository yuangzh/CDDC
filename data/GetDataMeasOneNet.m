function [y] = GetDataOneNet_y(G)
rand('seed',0); randn('seed',0);
[m,n] = size(G);
totnnz = round(0.1*n);
x = generate_x(n,totnnz);
b = G*x;
y = b + 0.1*randn(m,1);
y = max(0,y);



function x = generate_x(n,k)
p = randperm(n);
x = zeros(n,1);
x(p(1:k)) = randn(k,1);