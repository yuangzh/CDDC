function [x] = GetAvgCell(his)
% Example:
% Input:  his{1} = [1 3 4 5];
%         his{2} = [4 6 3 8 12 100];
%         his{3} = [5 7 0 5 1];
% Output:  [avg([1 4 5]) , avg([3 6 7]), avg([4 3 0]), avg([5 8 5]) ] = [3.33 5.33 2.33 6]
% Note that the length of x is min_i (his{i})

T = length(his);
N = inf;
for i=1:T
    N = min(N,length(his{i}));
end

for i = 1:N
    v = zeros(T,1);
    for j=1:T
        v(j) = his{j}(i);
    end
    x(i) = mean(v);
end