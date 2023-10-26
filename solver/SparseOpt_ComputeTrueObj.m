function [fobj] = SparseOpt_ComputeTrueObj(x,A,b,lambda,k)
fobj = 0.5*norm(A*x-b)^2 + lambda*norm(x,1) - lambda*topksum(x,k);