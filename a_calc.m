function a = a_calc(theta,fd,channels,pulse_length)
%A_CALC Summary of this function goes here
%   Detailed explanation goes here
wf=exp(-1i*(1:pulse_length)*2*pi*fd);
wt=exp(-1i*(0:(channels-1))*pi*sind(theta));
a=kron(wt,wf);
a=a/norm(a,2);
end

