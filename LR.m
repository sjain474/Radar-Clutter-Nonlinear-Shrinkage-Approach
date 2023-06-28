function RLS = LR(D,sigma_2,r)
%LR Summary of this function goes here
%   Detailed explanation goes here
[N,T]=size(D);
Lambda=ones(N,1);
D=D/sigma_2;
for k=1:r
    Lambda(k)=min([1 1/D(k)]);
end
RLS=sigma_2*(1./Lambda);
end

