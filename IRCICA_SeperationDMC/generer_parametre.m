function [A,C,D,F]=generer_parametre(p,alpha,N,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ici ce fonction sert a calculer les 4 parametre A C D F 
%formule 16 19 21 22 dans le papier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=2*(gamma(1-p)*cos(pi*p/2))/p;
F=2*gamma(1-p/alpha)/(gamma(1-p)*cos(pi*p/2));
eval(['C=(1/(2*pi)*quad(inline(''(abs(cos(x))).^(' num2str(alpha) ')'',''x''),-pi, pi))^(p/alpha);'])
%C=(quad(inline(['(abs(cos(x))).^(eval('alpha'))','x'),-pi, pi))^(p/alpha);
C=(1/(2*pi)*integral(@(x)(abs(cos(x))).^(1.83),-pi, pi)^(p/alpha));
n=(N-1)/(2*k)+1;
Hn=(pi*(2/pi)^(2*k*alpha)+(2*pi*k*alpha)/(2*k*alpha-1))*n^(2*k*alpha-1);
qkn=(1/2*(2/pi)^(2*k)+k/(2*k-1))*n^(2*k-1);
A=qkn/((Hn)^(1/alpha));
