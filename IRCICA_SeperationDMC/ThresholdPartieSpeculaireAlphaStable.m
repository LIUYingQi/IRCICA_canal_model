%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           ce ficher est dans le but de generer des simulation
%%%           pour tester le threhold N0 pour la simulation pour partie
%%%           speculaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      les formules indique dans ce programmme sont dans le papier
%%%      <<Statistical Channel Model Based on alpha-Stable Random Process
%%%      and Application tothe 60 GHz Ultra Wide Band Channel>>
%%%      Nourddine Azzaoui et Laurent Clavier, Menber IEEE, IEEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            load resultat dans simulation_IR.m
%        ce sont des parametre deja calcule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
Ntest=[1 6500]; %250 times * 26 point
alpha=1.83; %1.83 alpha 
p=alpha/2.4;
x=linspace(0,100,1601); 
preci=0.01;
risk=0.05;
pas=[0:0.5:100]; 
load('simulation_IR_result.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  generation des reponses impulsionelles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-definir N0 pour tester
N0_test =1:1:500;

% generer simulation
disp('generation des reponses impulsionelles avec N0 different')
omega=[1:0.0125:3];       
Hws=[];
f=1:0.00125:3;
dense_fh = mean(fh,1);

% pre-distribuer espace
amp=zeros(length(N0_test),1601);
AMP=zeros(length(N0_test),1601);

% pour trouver un N0 convient pour le threshold 
% pour chaque N0 ,on utilise la meme simulation et on fait 10 simulation
% choisir puissance moyen et delay spread comme parametre pour voir si N0 est bon ou pas.
% on essayer les N0 definit dans le debut

for j=1:10
    disp([ '---10  Generation--- : ' int2str(j)]);
    % processus de generer simulation avec differernt N0
    
    %%% mesure sur R en formule 25
    w=zeros(1,length(x)-1);
    for i= 1:(length(x)-1)
        w(i)=(x(i+1)-x(i))*max([dense_fh(i),dense_fh(i+1)]);
    end
    mE=sum(w);

    %%% C_alpha constant en formule 26
    C=(mE*(1-alpha)/(2^(alpha/2)*gamma(1+alpha/2)*gamma(2-alpha)*cos(pi*alpha/2)))^(1/alpha);

    %%% generer implusion reponse
    v=generer_v(x,N0,w);       % Instants
    h_const=generer_h_const(N0,alpha);   % Amplitude temp
    
    %%% pour chaque simulation on tester l'influence de N0
    for k=N0_test
        v_temp=v(1:k);
        amp(round(k/(N0_test(2)-N0_test(1)))+1,:)=generer_h(x,v_temp,h_const,C);
        u=(ones(length(f),1)*amp(round(k/(N0_test(2)-N0_test(1)))+1,:)).*exp(sqrt(-1)*f'*x);
        AMP(round(k/(N0_test(2)-N0_test(1)))+1,:)=C*sum(u.');
        AMP(round(k/(N0_test(2)-N0_test(1)))+1,:)=AMP(round(k/(N0_test(2)-N0_test(1)))+1,:)/sqrt(sum(abs(AMP(round(k/(N0_test(2)-N0_test(1)))+1,:)).^2));    
        amp(round(k/(N0_test(2)-N0_test(1)))+1,:)=abs(amp(round(k/(N0_test(2)-N0_test(1)))+1,:));
    end
    
    %%% sauvegarder amp et AMP pour cette simulation
    amp_nom = strcat('amp',int2str(j),'.mat');
    AMP_nom = strcat('AMP',int2str(j),'.mat');
    save(strcat('amp_temp',int2str(j),'.mat'),'amp');
    save(strcat('amp_freq',int2str(j),'.mat'),'AMP');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                calculer puissance moyenne et DS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

power = zeros(length(N0_test),10);
DS = zeros(length(N0_test),10);
power_result = zeros(1,length(N0_test));
DS_result = zeros(1,length(N0_test));
for k = N0_test
    tmp=floor(k/(N0_test(2)-N0_test(1)))+1;
    for j = 1:10
        amp_nom = strcat('amp_temp',int2str(j),'.mat');
        load(amp_nom);
        power(tmp,j) = mean(amp(tmp,:).^2);
        DS(tmp,j) = DelaySpread(amp(tmp,:),tsim,p);
    end
    power_result(tmp)=mean(power(tmp,:));
    DS_result(tmp)=mean(DS(tmp,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                comparaison des resultat obtenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(N0_test,power_result(2:end))   
title('profil puissance selon threshold');
xlabel('N0 pour threshold');
ylabel('puissace');
legend('puissace moyenne');
hold off

figure
hold on
plot(N0_test,DS_result(2:end))  
title('profil delay spread selon threshold');
xlabel('N0 pour threshold');
ylabel('DS');
legend('DS');
hold off