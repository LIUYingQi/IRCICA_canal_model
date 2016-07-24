%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           ce ficher est dans le but de generer fh indice 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      les formules indique dans ce programmme sont dans le papier
%%%      <<Statistical Channel Model Based on alpha-Stable Random Process
%%%      and Application tothe 60 GHz Ultra Wide Band Channel>>
%%%      Nourddine Azzaoui et laurent Clavier, Menber IEEE, IEEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

emplacement =10;   % changer ici pour emplacement
pos = 50:1:51;     % changer ici pour position (selon lamda length of wave we choose two )

napp=50;     % Number of impuls responses used for estimating mu.
Ntest=250;  % Number of available impulse responses 

alpha=1.83; % Chosen alpha - should be estimated...
p=alpha/2.4;
x=linspace(0,100,1601); % Time axis

preci=0.01; % If we setup the number of generated components to limit 
risl=0.05;  % the errorf

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Estimation of mu.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T indicates the impulse responses chosen for the estimation.   
T=pos; % We choose them no more randomly on the whole set.

% Chargement des donnes
fh=zeros(length(T),length(x));

for indi=1:1:length(T)
    %disp(['emplacement ' int2str(floor(T(indi)/250)+1) ', position ' int2str(mod(T(indi),250)+1)])
    disp('Enregistrement');
    disp([int2str(indi) ' /2']);
    disp(['emplacement 10 - position ' int2str(T(indi))]);
    H=[];
    while length(H)<=1
        % Enregistrement des fonctions de transfert
        % le while sert a verifier que le fichier existe bien, sinon le if
        % permet d'en changer.
        % [f,H]=Tcom_O9O9_loadH(mod(T(indi),250)+1,floor(T(indi)/250)+1);
        [f,H]=Tcom_O9O9_loadH(T(indi),10); % Je n'ai que l'emplacement 10.
        if length(H)<=1
            T=randi([1 Ntest],1,1);
        end
    end
    H=H/sqrt(sum(abs(H).^2));
    % f les frequences auquelles les fonctions de transfert sont observees.
    % H la vraie fonction de transfert mesure directement par le sondeur
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %				% Estimation de la densite spectrale %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre calcul des constantes de l'estimation
    disp('   ... Estimation');
    N=length(H);
    deltaf=2/(N-1);    % le pas d'echantionnage
    l=1;
    [A,C,D,F]=Tcom_O9O9_1_ctes(p,alpha,N,l);
    temp=0:0.01:1;
    dens0=(D*(A^alpha))/(F*C)*Tcom_O9O9_3_smooth(H,x,alpha,p,l,deltaf,N);
    fh(indi,:)=dens0;
end
mu=mean(fh,1);

save IndiceFHMesureChoisit.mat mu