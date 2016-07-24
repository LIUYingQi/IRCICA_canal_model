function [A,C,D,F]=Tcom_O9O9_1_ctes(p,alpha,N,k)

D=2*(gamma(1-p)*cos(pi*p/2))/p;
F=2*gamma(1-p/alpha)/(gamma(1-p)*cos(pi*p/2));
eval(['C=(1/(2*pi)*quad(inline(''(abs(cos(x))).^(' num2str(alpha) ')'',''x''),-pi, pi))^(p/alpha);'])
%C=(quad(inline(['(abs(cos(x))).^(eval('alpha'))','x'),-pi, pi))^(p/alpha);
%C=(1/(2*pi)*quad(inline('(abs(cos(x))).^(1.83)','x'),-pi, pi))^(p/alpha);

n=(N-1)/(2*k)+1;

Hn=(pi*(2/pi)^(2*k*alpha)+(2*pi*k*alpha)/(2*k*alpha-1))*n^(2*k*alpha-1);
qkn=(1/2*(2/pi)^(2*k)+k/(2*k-1))*n^(2*k-1);
A=qkn/((Hn)^(1/alpha));
