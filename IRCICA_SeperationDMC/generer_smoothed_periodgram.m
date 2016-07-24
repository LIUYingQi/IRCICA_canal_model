function f=generer_smoothed_periodgram(H,lambda,alpha,p,k,tau,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cette fonction est pour enfin calculer formule 24 
% Inputs
% H : la fonction de transfert (complexe, normalise, synchronise)
% lambda : points (instants) ou sont effectue les estimations
% alpha : exposant caracteistique
% p : coefficient (alpha/2.5 souvent)
% K : coefficient pour jackson's polynomial
% tau : coefficient pour Qn dans formule 13
% N : nombre de point de simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=-(N-1)/2:(N-1)/2;
q=generer_Qn(2*k,N);
hh=q(2*k,:)/q(2*k,N);   
% Un vecteur entre 0 et 1
% Ici hh va de 0 a 1 avec un pas de 1/N d'ou N=1601 est le nombre
% d'echantillons (frequentiels) dans la fonction de transfert
M=((N-1)/(2*k)+1)^(1/19);
a=lambda-1/M;
temp=0:0.01:1;

for i=1:length(lambda)
    u=a(i)+(2/M)*temp;
    h=(generer_IN(H,u,hh,N,alpha,m,tau)).^p;   % formule 17
    w=M*generer_epanechnikol_kernel(M*(lambda(i)-u));  % formule 23
    f(i)=mean(h.*w);
end
f=(2/M)*f.^(alpha/p);  % formule 24
