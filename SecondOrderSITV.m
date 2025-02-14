function [penalty] = SecondOrderSITV(f,epsilon)
B2f = SecondOrderDiff(f);
d = numel(f);
h = zeros(d,1);
for i = 1:d
    c = sqrt(B2f(i)^2+B2f(d+i)^2+B2f(2*d+i)^2+B2f(3*d+i)^2);
    if c > epsilon
        h(i) = c - epsilon/2;
    else
        h(i) = c^2/(2*epsilon);
    end
end    
penalty = sum(h);
end


