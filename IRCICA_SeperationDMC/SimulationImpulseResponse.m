%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           ce ficher est dans le but de generer 1000 simulation
%%%           pour tester si le modele est bon ou pas dans ce cas
%%%           dans la derniere partie plot deux dessin pour verifier le cas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      les formules indique dans ce programmme sont dans le papier
%%%      <<Statistical Channel Model Based on alpha-Stable Random Process
%%%      and Application tothe 60 GHz Ultra Wide Band Channel>>
%%%      Nourddine Azzaoui et Laurent Clavier, Menber IEEE, IEEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

disp('partie 1 definir des parametre')
Ntest=[1 6500]; %250 times * 26 point
napp=250;  %250 times pour chaque point
nsim=1000; %1000 simulation est fait pour le cas
alpha=1.83; %1.83 alpha 
p=alpha/2.4;
x=linspace(0,100,1601); 
preci=0.01;
risk=0.05;
pas=[0:0.5:100]; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%				Estimation de l'indice de stabilite alpha
%%%             (calculer l'indice fh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha=1.83;
% la m¨¦hode d'estimation de alpha n'est pas fiable
% p=alpha/2.5;
% reste a faire une mehode adaptative qui cherche le meilleur p.

disp('partie 2 Estimation de indice de stabilite alpha')
T=randi(Ntest,1,napp);
test=0;

% chercher 250 reponse impulse par hasard 
% c'est une facon pour generer random chiffre 
% 250 echantion pour estimation l'indice de alpha stable

while (test<=30)
    T=sort(T,'ascend');
    u=diff(T);
    z=find(u==0);
    if isempty(z)==0
        T(z)=randi(Ntest,1,length(z));
        test=test+1;
    else
        test=31;
    end
end

fh=zeros(length(T),length(x));
for indi=1:1:length(T)
    disp(['emplacement ' int2str(floor(T(indi)/250)+1) ', position ' int2str(mod(T(indi),250)+1)])
    disp('Chargement des donnes !!')
    H=[];
    while length(H)<=1
        [f,H]=loadH(mod(T(indi),250)+1,floor(T(indi)/250)+1);
        if length(H)<=1
            T=randi([1 Ntest],1,1);
        end
    end
    H=H/sqrt(sum(abs(H).^2));
    
% f les freuences auquelles les fonctions de transfert sont observes.
% H la vraie fonction de transfert mesure directement par le sondeur
% Estimation de la densite spectrale 
% pre-calcul des constantes de l'estimation

    N=length(H);       % le nombre d'echantionnage
    deltaf=2/(N-1);    % le pas d'ehantionnage (dans le papier on est 57-59 GHz)
    k=1;               % N = 2K(n-1)+1    
    [A,C,D,F]=generer_parametre(p,alpha,N,k);
    temp=0:0.01:1;
    density=(D*(A^alpha))/(F*C)*generer_smoothed_periodgram(H,x,alpha,p,k,deltaf,N);
    fh(indi,:)=density;   % fh matrice 250 * 1601
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  generation des reponses impulsionelles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('partie 3 generation des reponses impulsionelles')
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       delay spread
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('partie 4 delay spread')

%%%   delay spread pour simulation
p=2;
for i=1:nsim
    DS_sim(i)= DelaySpread(h(i,:),tsim,p);
end

%%%   delay spread pour mesure
Om_idx=1:500:1601;  
Hws=H(:,Om_idx);
Hw=[];
idx1=0;
DS_mesure=[];
for i=1:1:6499
    [f,H]=loadH(mod(i,250)+1,floor(i/250)+1);
    if isempty(H)
        idx1=idx1+1;
        H11(idx1,:)=H/sqrt(sum(abs(H).^2));
        Hw(idx1,:)=H(Om_idx)/sqrt(sum(abs(H).^2));
    end
    chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
    fichier=[chemin '\emplacement' int2str(floor(i/250)+1) '\C' int2str(mod(i,250)+1) 'temp.mat'];
    if exist(fichier)==0
        t=[];
        h=[];
        disp(['le fichier n''exsite pas - Emp ' int2str(floor(i/250)+1) ' Pos ' int2str(mod(i,250)+1)])
    else
        disp(['OK - Emp ' int2str(floor(i/250)+1) ' Pos ' int2str(mod(i,250)+1)])
        eval(['load ' fichier]);
        if (length(h_amp)==1601 && length(f)==1601)
            h=10.^(h_amp/20).*exp(sqrt(-1)*h_phase*pi/180); %.*exp(sqrt(-1)*(f)*2*pi*(distance(emplacement))/0.3/cos(20*pi/180));
            DS_mesure=[DS_mesure DelaySpread(h,t,p)];
        else
            disp(['Pb de taille - Emp ' int2str(floor(i/250)+1) ' Pos ' int2str(mod(i,250)+1)])
        end
    end
end

%%    plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fig 6 dans le papier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('partie 5 plot resultat')

%%%  voir le propriete de simulation et mesure pour delay spread pour voir
%%%  si la simulation et modele est correcte

%[f1,x1]=ksdensity(DS_sim,'function','cdf');
%[f2,x2]=ksdensity(DS_mesure,'function','cdf');

[f1,x1]=ksdensity(DS_sim);
[f2,x2]=ksdensity(DS_mesure);
