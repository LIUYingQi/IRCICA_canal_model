function z=generer_epanechnikol_kernel(x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEPAN       Multivariate Epanechikov Kernel Function.
% dans cette simulation on utilise KEPAN
% ce fonction est pour generer window Wn dans formule 23 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c(1)=2;
c(2)=pi;
if nargin==2,
  s=x.*x+y.*y;
else
  s=x.*x;
end;
z1=0.5*(nargin+2)*(1-s)/c(nargin);
z=z1.*(abs(s)<1);