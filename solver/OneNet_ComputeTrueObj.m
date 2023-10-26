function f = OneNet_ComputeTrueObj(x,G,y)
f = 0.5*norm(max(0,G*x) - max(0,y) )^2;
