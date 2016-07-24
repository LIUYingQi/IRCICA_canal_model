function f=Tcom_O9O9_3_smooth(H,lambda,alpha,p,k,tau,N)
% Inputs
% H : la fonction de transfert (complexe, normalisée, synchronisée)
% lambda : points (instants) où sont effectués les estimations
% alpha : exposant caractéristique
% p : coefficient (alpha/2.5 souvent)
% M : coefficient ???

m=-(N-1)/2:(N-1)/2;
q=Tcom_O9O9_2_ctes(2*k,N);
hh=q(2*k,:)/q(2*k,N);   % Un vecteur entre 0 et 1 (droite !)
  %%% Ici hh va de 0 à 1 avec un pas de 1/N où N=1601 est le nombre
  %%% d'achantillons (fréquentiels) dans la fonction de transfert
M=((N-1)/(2*k)+1)^(1/19);
a=lambda-1/M;
temp=0:0.01:1;

%dens0=(D*(A^alpha))/(F*C)*(2/M)*Tcom_O9O9_3_smooth(H,x,alpha,p,M,a,hh,m,deltaf,N,temp);

for i=1:length(lambda)
    u=a(i)+(2/M)*temp;
    h=(Tcom_O9O9_3_IN(H,u,hh,N,alpha,m,tau)).^p;
    w=M*Tcom_O9O9_3_Kernel(M*(lambda(i)-u));
    f(i)=mean(h.*w);
end
f=(2/M)*f.^(alpha/p);
