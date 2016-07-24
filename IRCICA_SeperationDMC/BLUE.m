clear
clc
close all
N = 5;
Mf = 1601;
bw = 2e9;
df = bw/(Mf-1);
dt = 1/bw;
delay = (0:Mf-1)*dt;
freq = (-(Mf-1)/2:(Mf-1)/2);
tau = 30e-9+60e-9*rand(N,1);
gamma = randn(N,1)+1i*randn(N,1);
gamma = gamma./max(abs(gamma));
B = exp(-2*pi*1i*freq'*df*tau')/sqrt(Mf);
X = B*gamma;
Br = sqrt(0.01/Mf)*(sqrt(2)/2)*(randn(Mf,1)+1i*randn(Mf,1));
H = X+Br;

plot(delay*1e9,20*log10(abs(ifft(H))))
gammas = B'*H;
db(abs(gamma-gammas))
(angle(gamma)-angle(gammas))*180/pi
Xs = B*gammas;
Hn = H - Xs;
hold on
%plot(delay*1e9,20*log10(abs(ifft(Hn))),'r')
plot(delay*1e9,20*log10(abs(ifft(Xs))),'k')
hold off

figure 
hold on
plot(abs(H));
plot(abs(Xs));
hold off
