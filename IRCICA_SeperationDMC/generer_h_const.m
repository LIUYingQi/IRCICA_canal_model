function h=generer_h_const(n,alp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generation une simulation de reponse impulsionnel
%%% parametre h_const en formule 28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp=-log(rand(1,n));
b=(cumsum(tmp)).^(-1/alp);
% the arrival times of poisson process
% distributed as gamma of parameter i
h=b.*(randn(1,n)+sqrt(-1)*randn(1,n)); 
%ajouter unit circle
