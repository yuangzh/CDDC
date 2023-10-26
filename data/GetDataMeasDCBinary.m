function [y] = GetDataMeasDCBinary(G)
rand('seed',0);
randn('seed',0);
[m,n] = size(G);
x = sign(randn(n,1));
b = G*x;
y = b + randn(m,1);
