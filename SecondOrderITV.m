function [penalty] = SecondOrderITV(f)
Bf = SecondOrderDiff(f);
d = numel(f);
%N = sqrt(d);
p = zeros(d,1);
for i = 1:d
    p(i) = sqrt(Bf(i)^2+Bf(d+i)^2+Bf(2*d+i)^2+Bf(3*d+i)^2);
end

%------------------------------------------------------------------
%Full Image Fisrt-Order ITV Penalty
%------------------------------------------------------------------
penalty = sum(p);
%------------------------------------------------------------------
end