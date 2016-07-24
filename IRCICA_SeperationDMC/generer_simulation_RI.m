function [amp,pas,AMP]=generer_simulation_RI(x,alpha,fh,N0,pas,f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cette fonction est pour generer simulation pour reponse impulsionelle
%%% inputs :
%%% x: points d'echnation
%%% alpha : parametre alpha stable
%%% fh : density qu'on a estime
%%% N : nombre de point
%%% pas : pas 
%%% f : frequence dans le periodgram on va generer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% mesure sur R en formule 25
w=zeros(1,length(x)-1);
for i= 1:(length(x)-1)
    w(i)=(x(i+1)-x(i))*max([fh(i),fh(i+1)]);
end
mE=sum(w);

%%% C_alpha constant en formule 26
C=(mE*(1-alpha)/(2^(alpha/2)*gamma(1+alpha/2)*gamma(2-alpha)*cos(pi*alpha/2)))^(1/alpha);

%%% generer implusion reponse
v=generer_v(x,N0,w);       % Instants
h_const=generer_h_const(N0,alpha);   % Amplitude temp
amp=generer_h(pas,v,h_const,C);

u=(ones(length(f),1)*amp).*exp(sqrt(-1)*f'*pas);
AMP=C*sum(u.');
AMP=AMP/sqrt(sum(abs(AMP).^2));    % Amplitude freq

