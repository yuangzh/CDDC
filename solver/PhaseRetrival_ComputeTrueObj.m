function [fobj] = PhaseRetrival_ComputeTrueObj(v,Q)
% min_{v} 0.5 v'Qv, s.t. v \in {-1,+1}^n
v=ProjBinary(v);
fobj = 0.5*v'*Q*v;

