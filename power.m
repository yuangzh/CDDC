n = 10;
A = randn(n);
x = randn(n,1);
for iter = 1:100
    y = A*x;
    x = y / norm(y); 
    norm(A*x / (x'*A*x) - x)
end

x

