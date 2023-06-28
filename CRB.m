clear all
close all
clc
%% Challenge DataSet
load('dataset_san_diego_scenario2_CPI11.mat')
Hc=reshape(true_clutter_channel_cube_san_diego_s2(1:8,:,:),64*8,[]);
x=waveform_san_diego_s2;
X=convmtx(x,size(Hc,2));
for j=1:size(Hc,1)
    HCs(:,j)=X*transpose(Hc(j,:));%Clutter Matrix
end
channels=8;
pulse_length=64;
Hcs=transpose(HCs);
[N1,T]=size(Hcs);
var=5e-14;%Noise Variance
S1=(1/T)*(Hcs*Hcs')+var*eye(512);
[Us1,Ds1,Vs1]=svd(S1);
invS1=pinv(S1);
b=floor(T/N1);
theta=30;
fd=0.2;
%%
MonteCarlo=1024;
for c=1:MonteCarlo
for j=1:b
Datapoints=j*N1;
inData=mvnrnd(zeros(1,N1),S1,Datapoints);
x=(randn(N1,Datapoints)+1i*randn(N1,Datapoints))/sqrt(2);
RSR=(x*(x'))/Datapoints;
[U,D,V]=svd(RSR);
LR_0=sum(log(diag(D)))+N1-trace(D);
gamma=N1/Datapoints;
S=(1/Datapoints)*inData'*(inData);
invS=pinv(S);
[Us,Ds,Vs]=svd(S);
ds=diag(Ds);
[eigenvals_stein, sig_stein, timeElapsedDon] = ...
    optimal_shrinkage(diag(Ds),gamma,'Stein');
[eigenRCML,sigma_rcml] = RCML(Ds,LR_0);
sigmahat_don_stein=Us*diag(eigenvals_stein)*Vs';
inv_sigmahat_don_stein=Us*diag(1./eigenvals_stein)*Vs';
inv_RCML=Us*diag(1./eigenRCML)*Vs';
alpha_max=ds(150);
alpha_min=ds(100);
inv_DLmax=Us*diag(1./(ds+alpha_max*ones(size(ds))))*Vs;
inv_DLmin=Us*diag(1./(ds+alpha_min*ones(size(ds))))*Vs;
for i=1:length(theta)
    for k=1:length(fd)
a = a_calc(theta(i),fd(k),channels,pulse_length);
SINR_est(i,k,c,j)=1/abs(a*inv_sigmahat_don_stein*a');
SINR_RCML(i,k,c,j)=1/abs(a*inv_RCML*a');
SINR_DLmax(i,k,c,j)=1/abs(a*inv_DLmax*a');
SINR_DLmin(i,k,c,j)=1/abs(a*inv_DLmin*a');
SINR_true(i,k)=1/abs(a*invS1*a');
    end
end
end
end
SINR_data=squeeze(mean(SINR_est,[1,2,3]));
SINR_Azimuth=squeeze(mean(SINR_est(:,:,:,2),[2,3]));
SINR_Doppler=squeeze(mean(SINR_est(:,:,:,2),[1,3]));
%%
SINR_data_RCML=squeeze(mean(SINR_RCML,[1,2,3]));
%%
SINR_data_DLmax=squeeze(mean(SINR_DLmax,[1,2,3]));
%%
SINR_data_DLmin=squeeze(mean(SINR_DLmin,[1,2,3]));
%%
SINR_data_true=squeeze(mean(SINR_true,"all"));
%%
figure(1)
plot(N1*(1:b),(1/var)*SINR_data,LineWidth=2,LineStyle="-",Marker="square",MarkerSize=8);
hold on;
plot(N1*(1:b),(1/var)*SINR_data_RCML,LineWidth=2,LineStyle="-.",Marker="*",MarkerSize=8);
hold on;
plot(N1*(1:b),(1/var)*SINR_data_true*ones(b,1),LineWidth=2,LineStyle="--",Marker=">",MarkerSize=8);
hold on;
plot(N1*(1:b),(1/var)*SINR_data_DLmax,LineWidth=2);
hold on;
plot(N1*(1:b),(1/var)*SINR_data_DLmin,LineWidth=2);
xlabel('Datapoints',FontSize=15,Interpreter='latex');
ylabel('Normalized Error Variance',Fontsize=15,Interpreter='latex');
grid on;
legend('Proposed','RCML-EL','True',"$\textnormal{DL}_{\textnormal{max}}$","$\textnormal{DL}_{\textnormal{min}}$",'Location','southeast','Interpreter','latex',Fontsize=14);
xlim([N1,N1*b]);
ylim([1,inf]);
set(gca,'FontSize',12);