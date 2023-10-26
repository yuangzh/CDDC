function [fobj] = BinaryOpt_ComputeTrueObj(x,G,y,lambda,n)
x = max(-1,min(1,x));
fobj = 0.5* norm(G*x-y)^2  + lambda*(sqrt(n) - norm(x));