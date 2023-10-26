function [x] = prox_l1(a,lambda)
% solving the following OP:
% min_{x} 0.5||x-a||^2 + lambda * sum(abs(x))
x = sign(a).*max(0,abs(a)-lambda);
 
