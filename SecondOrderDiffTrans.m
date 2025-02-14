function [BTb] = SecondOrderDiffTrans(b)
%b has to be a vector
d = numel(b)/4;
N = sqrt(d);
b1 = reshape(b(1:d),[N N]);
b2 = reshape(b(d+1:2*d),[N N]);
b3 = reshape(b(2*d+1:3*d),[N N]);
b4 = reshape(b(3*d+1:4*d),[N N]);
BTb1 = zeros(d,1);
BTb2 = zeros(d,1);
BTb3 = zeros(d,1);
BTb4 = zeros(d,1);
%% Calculate D_xx^T*b1(=D_xx*b1)
for i = 1:N
    BTb1((i-1)*N+1) = b1(2,i)-b1(1,i);
    BTb1(i*N) = b1(N-1,i)-b1(N,i);
    for j = 2:N-1
        BTb1((i-1)*N+j) = b1(j-1,i)+b1(j+1,i)-2*b1(j,i);
    end
end

%% Calculate D_xy^T*b2(=D_yx*b2)
for i=2:N
    BTb2((i-1)*N+1) = b2(2,i)-b2(2,i-1);
    BTb2(i*N) = b2(N,i-1)-b2(N,i);
    for j=2:N-1
        BTb2((i-1)*N+j) = (b2(j,i-1)-b2(j,i))-(b2(j+1,i-1)-b2(j+1,i));
    end
end


%% Calculate D_yx^T*b3(=D_xy*b2)
for j=2:N
    BTb3(j) = b3(j,2)-b3(j-1,2);
    BTb3(d-N+j) = b3(j-1,N)-b3(j,N);
end
for i = 2:N-1
    for j = 2:N
        BTb3((i-1)*N+j) = (b3(j-1,i)-b3(j,i))-(b3(j-1,i+1)-b3(j,i+1));
    end
end

%% Calculate D_yy^T*b4(=D_yy*b4)
for j=1:N
    BTb4(j) = b4(j,2)-b4(j,1);
    BTb4(d-N+j) = b4(j,N-1)-b4(j,N);
end
for i=2:N-1
    for j=1:N
        BTb4((i-1)*N+j) = b4(j,i-1)+b4(j,i+1)-2*b4(j,i);
    end
end

BTb = BTb1+BTb2+BTb3+BTb4;

end