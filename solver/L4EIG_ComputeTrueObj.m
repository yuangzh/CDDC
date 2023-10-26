function fobj = L4EIG_ComputeTrueObj(x,A,beta)
x = x ./ norm(x);
fobj = x'*A*x + beta*norm(x,4)^4;