function [BTb] = FirstOrderDiffTrans(b)
d = numel(b)/2;
N = sqrt(d);
ub = reshape(b(1:d),[N N]);
lb = reshape(b(d+1:2*d),[N N]);
DTub = zeros(d,1);
Ilb = zeros(d,1);
%BTb = DTub + Ilb

%compute DTub
for i = 1:N
    DTub((i-1)*N+1) = -ub(2,i);
    DTub(i*N) = ub(N,i);
end
for i = 1:N
    for j = 2:N-1
        DTub((i-1)*N+j) = ub(j,i) - ub(j+1,i);
    end
end

%compute Ilb
for j = 1:N
    Ilb(j) = -lb(j,2);
    Ilb((N-1)*N+j) = lb(j,N);
end
for i = 2:N-1
    for j = 1:N
        Ilb((i-1)*N+j) = lb(j,i) - lb(j,i+1);
    end
end

BTb = DTub + Ilb;
end
    