function [eigenRCML,sigma_rcml] = RCML(Ds,LR_0)
%RCML Summary of this function goes here
%   Detailed explanation goes here
r_0=54;
[N,~]=size(Ds);
ds=diag(Ds);
sigma_2=zeros(100,1);
l=1;
sigma_2(l)=(1/r_0)*sum(diag(Ds(N-r_0:N,N-r_0:N)));
rank=r_0;
err=1e-6;
while(err>1e-19)
while(rank<N && rank>0)
    dd=sort(LR(diag(Ds),sigma_2(l),rank),'descend');
    ddp1=sort(LR(diag(Ds),sigma_2(l),rank+1),'descend');
    ddm1=sort(LR(diag(Ds),sigma_2(l),rank-1),'descend');
    LRp1=(sum(log(ds./ddp1))+N-sum(ds./ddp1)-LR_0)^2;
    LRm1=(sum(log(ds./ddm1))+N-sum(ds./ddm1)-LR_0)^2;
    LR1=(sum(log(ds./dd))+N-sum(ds./dd)-LR_0)^2;
    if(LR1<=LRp1 && LR1<=LRm1)
        break;
    end
    if(LR1<=LRp1 && LR1>=LRm1)
            rank=rank-1;
    end
    if(LR1>=LRp1 && LR1<=LRm1)
            rank=rank+1;
    end
end
r_0=rank;
sigma_2(l+1)=1/(N-r_0)*sum(ds(r_0+1:N));
err=abs(sigma_2(l+1)-sigma_2(l));
l=l+1;
end
eigenRCML=dd;
sigma_rcml=sigma_2(100);
end

