function Seuil=generer_threshold(t,preci,risk,alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The treshold of the infinite sum giving an approximation with 
% precision preci and risk probability risk.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Seuil=floor((length(t)/2)^(1/alpha)/preci/risk)+1;