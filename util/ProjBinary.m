function [v] = ProjBinary(v)
% min_v ||v - a||_2^2, s.t. v \in {-1,+1}^n
absv = abs(v);
v = v./absv;
v(absv==0) = 1;