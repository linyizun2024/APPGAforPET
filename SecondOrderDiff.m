function [Bf] = SecondOrderDiff(f)
%f has to be a vector
d = numel(f);
N = sqrt(d);
f = reshape(f,[N N]);
Bf = zeros(4*d,1);
%% Calculate D_xx: [I\Kro((-D^T)D)](f)
for i = 1:N
    Bf((i-1)*N+1) = f(2,i)-f(1,i);
    Bf(i*N) = f(N-1,i)-f(N,i);
    for j = 2:N-1
        Bf((i-1)*N+j) = f(j-1,i)+f(j+1,i)-2*f(j,i);
    end
end

%% Calculate D_xy: [(-D^T)\KroD](f)
for j=2:N
    Bf(d+j) = f(j,2)-f(j-1,2);
    Bf(2*d-N+j) = f(j-1,N)-f(j,N);
end
for i = 2:N-1
    for j = 2:N
        Bf(d+(i-1)*N+j) = (f(j-1,i)-f(j,i))-(f(j-1,i+1)-f(j,i+1));
    end
end

%% Calculate D_yx: [D\Kro(-D^T)](f)
for i=2:N
    Bf(2*d+(i-1)*N+1) = f(2,i)-f(2,i-1);
    Bf(2*d+i*N) = f(N,i-1)-f(N,i);
    for j=2:N-1
        Bf(2*d+(i-1)*N+j) = (f(j,i-1)-f(j,i))-(f(j+1,i-1)-f(j+1,i));
    end
end

%% Calculate D_yy: [((-D^T)D)\KroI](f)
for j=1:N
    Bf(3*d+j) = f(j,2)-f(j,1);
    Bf(4*d-N+j) = f(j,N-1)-f(j,N);
end
for i=2:N-1
    for j=1:N
        Bf(3*d+(i-1)*N+j) = f(j,i-1)+f(j,i+1)-2*f(j,i);
    end
end
    
end