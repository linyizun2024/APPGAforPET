function [disk] = FieldofView(nR)
nn = (nR-1)/2;
x = linspace(-nn,nn,nR).^2;
disk = x(ones(nR,1),:) + x(ones(1,nR),:)' <= (nR/2)^2;
disk = double(disk);
end