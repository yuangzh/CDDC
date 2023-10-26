function [t_best] = ncvx_prox_poly4(c1,c2,c3,c4,lb,ub)
% min c1*t + c2*t*t + c3*t*t*t + c4*t*t*t*t
% grad = c1 + 2*c2*t + 3*c3*t*t + 4*c4*t*t*t = 0
xs = roots([4*c4,3*c3,2*c2,c1]);
xs = real(xs);
xs = min(ub,xs);
xs = max(lb,xs);
xs = [xs;lb;ub];
for i=1:length(xs)
    HandleObj = @(t)c1*t + c2*t*t + c3*t*t*t + c4*t*t*t*t;
    fs(i) = HandleObj(xs(i));
end
[~,index]=min(fs);
t_best = xs(index);