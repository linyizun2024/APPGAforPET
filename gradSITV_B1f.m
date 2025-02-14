function gradg_atB1f = gradSITV_B1f(B1f,epsilon)
d = numel(B1f)/2;
gradg_atB1f = zeros(2*d,1);
for i = 1:d
    c = max(sqrt(B1f(i)^2+B1f(d+i)^2),epsilon);
    gradg_atB1f(i) = B1f(i)/c;
    gradg_atB1f(i+d) = B1f(i+d)/c;
end
end
