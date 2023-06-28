function [D_k,U_k,h_k] = rank_calc(R,sigma2,gamma)
%RANK_CALC Summary of this function goes here
%   Detailed explanation goes here
[U,D,V]=svd(R);
d=diag(D);
k=sum((d>sigma2*(1+sqrt(gamma))));
D_k=d(1:k)-sigma2*ones(size(d(1:k)));
h_k=(D_k.^2-sigma2^2*gamma*ones(size(d(1:k))))./(D_k.*(D_k-sigma2*gamma*ones(size(d(1:k)))));
U_k=zeros(size(R,1),k);
for j=1:k
    U_k(:,j)=U(:,j)/norm(U(:,j),2);
end
end
