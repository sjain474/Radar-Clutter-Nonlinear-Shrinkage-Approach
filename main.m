clear all
close all
%% Challenge DataSet
load('dataset_san_diego_scenario2_CPI11.mat')
Hc=reshape(true_clutter_channel_cube_san_diego_s2(1:8,:,:),64*8,[]);
x=waveform_san_diego_s2;
X=convmtx(x,size(Hc,2));
for j=1:size(Hc,1)
    HCs(:,j)=X*transpose(Hc(j,:));
end
Hcs=transpose(HCs);
var=5e-14;%Noise Variance
n=(randn(size(Hcs))+1i*randn(size(Hcs)))*sqrt(var/2);
y_l=Hcs+n;
%% For Challenge DataSet
y=y_l(:,:);
yc=Hcs;
yn=n;
%%
[N1,T]=size(y);
Ta=size(yc,2);
Tn=size(yn,2);
S1=(1/T)*y*(y');
gamma=0.5;
Datapoints=N1/gamma;
ys=y(:,1:Datapoints);
inData=ys;
[N,Ts]=size(inData);
S=(1/Ts)*inData*(inData');
[Us,Ds,Vs]=svd(S);
ds=diag(Ds);
[eigenvals_stein, sig_stein, timeElapsedDon] = optimal_shrinkage(diag(Ds),gamma,'Stein');
sigmahat_don_stein=Us*diag(eigenvals_stein)*Vs';
inv_sigmahat_don_stein=Us*diag(1./eigenvals_stein)*Vs';
%%
[D_k,U_k,h_k] = rank_calc(S1,var,gamma);
channels=8;
pulse_length=64;
theta=30;
fd=0.2;
a = a_calc(theta,fd,channels,pulse_length);
SINR=1/abs(a*inv_sigmahat_don_stein*a');
for i=1:length(theta)
    for j=1:length(fd)
[v_0(i,j),v_1(i,j)] = v_0v_1(U_k,D_k,theta(i),fd(j),channels,pulse_length,h_k,var);
    end
end
%%
for k=1:5
pfa(k)=1-10^(-k);    
x=chi2inv(pfa(k),2);
snr_dB=-10:2:30;
snr=10.^(snr_dB/10);
for l=1:length(snr)
for i=1:length(theta)
    for j=1:length(fd)
delta=snr(l)/v_1(i,j);
alpha2=var*snr(l);
sinr(l)=10*log10(alpha2/SINR);
P(i,j)=ncx2cdf(x,2,delta);
    end
end
Pd(l,k)=1-sum(P,'all')/(length(fd)*length(theta));
legendCell{k}="pfa="+num2str(10^(-k));
end
end
%%
plot(sinr,Pd);
ylabel('$P_d$','FontSize',15,'Interpreter','latex');
xlabel('SNR (dB)','FontSize',15,'Interpreter','latex');
legend(legendCell,'location','Southeast','fontsize',15);
xlim([sinr(1),sinr(length(sinr))]);
grid on