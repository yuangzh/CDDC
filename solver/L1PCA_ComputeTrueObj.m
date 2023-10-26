function f = L1PCA_ComputeTrueObj(x,G)
f = 0.5*x'*x - norm(G*x,1) ;
