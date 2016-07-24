%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ce programme g¨¦n¨¨re une DMC par deux m¨¦thodes
%%% 1ere m¨¦thode: D¨¦composition de Cholesky
%%% 2¨¨me m¨¦thode: utilise structure de Toeplitz
%%% Version Debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars;
% close all;
clc;

%% PARAMETRES SIMULATIONS
Mf         = 1601;
N          = 1000;

tol_a0      = 1e-6;
tol_a1      = 1e-6;
tol_bd      = 1e-6;
tol_td      = 1e-6;

%% PARAMETRES DMC
a0_dB     = -87;
a1_dB     = -54;
a0        = 10^(a0_dB/10);
a1        = 10^(a1_dB/10);
beta_d    = 0.005;
tau_d     = 0.1;

%% CONSTRUCTION MODELE
ka       = transpose((a1/Mf)*(exp(-1i*2*pi.*(0:Mf-1)*(tau_d))./(beta_d+1i*2*pi.*(0:Mf-1)/Mf)));
ka(1)    = ka(1) + a0;
Rf_dmc   = toeplitz(ka,ka');

%% TIRAGE DMC METHODE DE CHOLESKY
z        = sqrt(2)/2*(randn(Mf,N)+1i*randn(Mf,N));
L        = chol(Rf_dmc,'lower');
Xf_I     = L*z;

% F        = zeros(Mf,Mf);
% w        = exp(-1i*2*pi/Mf);
% for j = 0 : Mf-1
%     for i = 0 : Mf-1
%         F(i+1,j+1) = w^(i*j);
%     end
% end
% F   = F./sqrt(Mf);
% clear w
% Xt_I       = (1/N)*sum((F'*Xf_I).*conj(F'*Xf_I),2);

Xt_I = mean((abs(ifft(Xf_I)).^2)*Mf,2);

figure
plot(10*log10(abs(Xt_I)),'r')
