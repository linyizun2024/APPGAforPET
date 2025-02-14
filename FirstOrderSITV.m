function [penalty] = FirstOrderSITV(f,epsilon)
B1f = FirstOrderDiff(f);
d = numel(f);
h = zeros(d,1);
for i = 1:d
    c = sqrt(B1f(i)^2+B1f(d+i)^2);
    if c > epsilon
        h(i) = c - epsilon/2;
    else
        h(i) = c^2/(2*epsilon);
    end
end    
penalty = sum(h);
end

