function [v_0,v_1] = v_0v_1(U_k,D_k,theta,fd,channels,pulse_length,h_k,sigma2)
%V_0V_1 Summary of this function goes here
%   U_k: A matrix of size pxk
%   D_k: An array that contains the eigenvalues lambda-sigma^2
%   theta: target azimuth
%   fd: Target Doppler
%   channels: No.of Array elements
%   pulse_length: Length of the pulses
%   h_k: The cosine between the leading eigenvectors
%   sigma2:Noise power
a = a_calc(theta,fd,channels,pulse_length);
[N,k]=size(U_k);
Pi=eye(N)-U_k*U_k';
num=zeros(k,1);
denom=zeros(k,1);
for j=1:k
    num(j)=(D_k(j)/sigma2)*((1-h_k(j))^2)*abs(a*U_k(:,j))^2;
    denom(j)=((1-h_k(j))^2)*abs(a*U_k(:,j))^2;
end
v_0=1+(sum(num))/(norm(Pi*a',2)^2+sum(denom));
v_1=v_0/(norm(Pi*a',2)^2+sum(denom));
end
