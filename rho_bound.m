clear all
close all
%% Challenge DataSet
load('dataset_san_diego_scenario2_CPI11.mat')
Ht=reshape(true_target_channel_cube_san_diego_s2(1:8,:,:),64*8,[]);%Clutter Contaminated With Target
x=waveform_san_diego_s2;
X=convmtx(x,size(Hc,2));
for j=1:size(Hc,1)
    HCs(:,j)=X*(transpose((Ht(j,:))+Ht(j,:)));
end
Hcs=transpose(HCs);
%%
M=1024;
[N1,T]=size(Hcs);
b=floor(T/N1);
%%
fd=-0.5:0.05:0.5;
theta=-180:5:180;
rho_don_Stein=zeros(length(fd),length(theta),b,M);
rho_rie=zeros(length(fd),length(theta),b,M);
for m=1:M
var=5e-14;%Noise Variance
n=(randn(size(Hcs))+1i*randn(size(Hcs)))*sqrt(var/2);
y_l=Hcs+n;
y=y_l(:,:);
yc=Hcs;
yn=n;
[N1,T]=size(y);
Ta=size(yc,2);
Tn=size(yn,2);
S1=(1/T)*y*(y');
invS1=pinv(S1);
A=(1/Ta)*yc*(yc');
Ns=(1/Tn)*yn*(yn');
[Ua,Da,Va]=svd(A);
[U,D,V]=svd(S1);
[Un,Dn,Vn]=svd(Ns);
d=diag(D);
dn=diag(Dn);
da=diag(Da);
eigenvals(:,m)=d/var;
%%
for j=1:b
Datapoints=j*N1;
gamma=N1/Datapoints;
I = (eigenvals(:,m) > 1+sqrt(gamma));
eigenvals(~I,m) = 1;
ys=y(:,1:Datapoints);
inData=ys;
[N,Ts]=size(inData);
x=(randn(N,Ts)+1i*randn(N,Ts))/sqrt(2);
RSR=(x*(x'))/Ts;
[U,D,V]=svd(RSR);
LR_0=sum(log(diag(D)))+N-trace(D);
S=(1/Ts)*inData*(inData');
[Us,Ds,Vs]=svd(S);
ds=diag(Ds);
%% RCML EL
tic
r_0=54;
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
l=l+1
end
sigmahat_rie=Us*diag(dd)*Vs';
timeElapsedRIE(j,m)=toc;
inv_sigmahat_rie=Us*diag(1./dd)*Vs';
l=1;
%% Proposed Shrinkage
[eigenvals_stein, sig_stein, timeElapsedDon ] = optimal_shrinkage(diag(Ds),gamma,'Stein');
eta(:,j,m)=Stein(eigenvals(:,m),gamma);
sigmahat_don_stein=Us*diag(sqrt(eigenvals_stein))*Vs';
inv_sigmahat_don_stein=Us*diag(1./sqrt(eigenvals_stein))*Vs';
k_delta(j,m)=cond(squeeze(eta(:,j,m)),squeeze(eigenvals(:,m)),gamma);
rho_min(j,m)=10*log10(4/(k_delta(j,m)+(1/k_delta(j,m))+2));
for k=1:length(fd)
    for l=1:length(theta)
wf=exp(-1i*(1:64)*2*pi*fd(k));
        wt=exp(-1i*(0:7)*pi*sind(theta(l)));
w=kron(wf,wt);
rho_rie(k,l,j,m)=abs((w*inv_sigmahat_rie*w')^2/((w*invS1*w')*(w*inv_sigmahat_rie*S1*inv_sigmahat_rie*w')));
rho_don_Stein(k,l,j,m)=abs((w*inv_sigmahat_don_stein*w')^2/((w*invS1*w')*(w*inv_sigmahat_don_stein*S1*inv_sigmahat_don_stein*w')));
    end
end
end
end
%%
rho_Rie=squeeze(10*log10(mean(mean(mean(rho_rie,4),1),2)));
rho_Stein=squeeze(10*log10(mean(mean(mean(rho_don_Stein,4),1),2)));
rho_m=mean(rho_min,2);
%% Plots
figure();
plot(N1*(1:(b-1)),rho_m(1:(b-1)));
hold on
plot(N1*(1:(b-1)),rho_Stein(1:(b-1)),'x-');
hold on;
plot(N1*(1:(b-1)),rho_Rie(1:(b-1)),'x-');
hold on;
plot(N1*(1:(b-1)),zeros(size((1:(b-1)))));
xlabel('Datapoints','fontsize',14);
ylabel('\rho(dB)','fontsize',14);
legend('Lower Bound','Stein Loss','RCML-EL','Upper Bound','fontsize',12,'location','Southeast');
grid on;
ylim([-3.5,0.1]);
xlim([N1,(b-1)*N1]);
%% Functions
function k_delta=cond(eta,l,gamma)
n=length(eta);
c_2=cosine(l,gamma);
s_2=ones(n,1)-c_2;
[v_plus, v_minus]=quadratic(eta,l,c_2,s_2);
k_delta=max([1;v_plus])/min([1;v_minus]);
end
function D=determinant(eta,l)
D=eta./l;
end
function T=Trace(eta,l,c,s)
T=(((s+c.*eta)./l)+c+s.*eta);
end
function [v_plus, v_minus]=quadratic(eta,l,c,s)
T=Trace(eta,l,c,s);D=determinant(eta,l);
v_plus=T/2+sqrt((T/2).^2-D);
v_minus=T/2-sqrt((T/2).^2-D);
end
function c_2=cosine(l,gamma)
n=length(l);
c_2=zeros(n,1);
I=(l>1+sqrt(gamma));
c_2(I)=(ones(length(l(I)),1)-(gamma*ones(length(l(I)),1))./(l(I)-ones(length(l(I)),1)))...
    ./(ones(length(l(I)),1)+(gamma*ones(length(l(I)),1))./(l(I)-1).^2);
end
function eta_stein=Stein(l,gamma)
c_2=cosine(l,gamma);
s_2=ones(length(c_2),1)-c_2;
I=(l>1+sqrt(gamma));
eta_stein=ones(length(l),1);
eta_stein(I)=l(I)./(c_2(I)+s_2(I).*l(I));
end