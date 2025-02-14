function [penalty] = FirstOrderITV(f)
Bf = FirstOrderDiff(f);
d = numel(f);

p = zeros(d,1);
for i = 1:d
    p(i) = sqrt(Bf(i)^2+Bf(d+i)^2);
end
%------------------------------------------------------------------
penalty = sum(p);
%------------------------------------------------------------------
end