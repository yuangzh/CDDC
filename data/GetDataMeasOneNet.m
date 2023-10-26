function [y] = GetDataMeasOneNet(G)
[m,n] = size(G);
x = randn(n,1);
b = G*x;
y = b + 0.1*norm(b)*randn(m,1);
y = max(0,y);






