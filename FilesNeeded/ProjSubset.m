function [Subset] = ProjSubset(SubsetNum,nPhi,nLOR,nSubsets)

m=nPhi/nSubsets;
Sub=zeros(m,1);
Subset=zeros(m*nLOR,1);
for i=1:m
    Sub(i)=(i-1)*nSubsets+SubsetNum;
end
for i=1:m
    for j=1:nLOR
        Subset((i-1)*nLOR+j)=(Sub(i)-1)*nLOR+j;
    end
end
