function [Bf] = FirstOrderDiff(f)
N = sqrt(numel(f));
f=reshape(f,[N N]); 
Bf = zeros(2*N*N,1);
for i = 1:N
    for j = 2:N
        Bf((i-1)*N+j) = f(j,i)-f(j-1,i);
    end
end
for i = 2:N
    for j = 1:N
        Bf((N+i-1)*N+j) = f(j,i)-f(j,i-1);
    end
end
end
    
