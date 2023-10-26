function [Q,normQ] = PhaseRetrivalComputeQ(G,y)
% W = G*pinv(G,1e-7)*diag(y)-diag(y);
W = G*pinv(full(G))*diag(y)-diag(y);
Q = W'*W ;
normQ = SpectralNorm(W) ;

