function s=generer_IN(X,lambda,h,N,alpha,m,tau)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ce fonction est pour generer formule 18 dans le papier
%  c'est pour en fin avoir expression IN(lamda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e=exp(sqrt(-1)*(((tau*m))'*(2*pi*lambda)));
s=(abs((tau^(1/alpha))*real(((h.*X)*e))));