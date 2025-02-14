function gradg_atB2f = gradSITV_B2f(B2f,epsilon)
d = numel(B2f)/4;
gradg_atB2f = zeros(4*d,1);
for i = 1:d
    c = max(sqrt(B2f(i)^2+B2f(i+d)^2+B2f(i+2*d)^2+B2f(i+3*d)^2),epsilon);
    gradg_atB2f(i) = B2f(i)/c;
    gradg_atB2f(i+d) = B2f(i+d)/c;
    gradg_atB2f(i+2*d) = B2f(i+2*d)/c;
    gradg_atB2f(i+3*d) = B2f(i+3*d)/c;
end
end

