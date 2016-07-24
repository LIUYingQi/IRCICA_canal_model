%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de comparer simulation
%%% et mesure deja fait.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  load fichier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
load ('simulation_IR_result.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  generation des simulation reponses impulsionelles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('partie generation simulation des reponses impulsionelles')
omega=[1:0.0125:3];       
Hws=[];
f=1:0.00125:3;
h=zeros(nsim,length(x));
H=zeros(nsim,length(f));
N0=3*generer_threshold(x,preci,risk,alpha)+10000; 

for j=1:nsim
    disp([' --- ' int2str(nsim) ' ---  Generation - ' int2str(j)])
    [h(j,:),tsim,H(j,:)]=generer_simulation_RI(x,alpha,mean(fh,1),N0,x,f);
end
h_simulation=h;      %    h_simulation reponse temporelle simulation   1000 simulation avec chaque simulation 1601 points
H_simulation=H;      %    H_simulation reponse spectralle simulation   1000 simulation avec chaque simulation 1601 points
