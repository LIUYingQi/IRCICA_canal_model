%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           ce ficher est dans le but de generer des simulation
%%%           pour voir l'influence de threshold N0 pour partie speculaire
%%%           et partie dense
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

emplacement =10;
pos =50:1:51;

napp=30;     % Number of impuls responses used for estimating mu.
Ntest=250;  % Number of available impulse responses 

alpha=1.83; % Chosen alpha - should be estimated...
p=alpha/2.4;
x=linspace(0,100,1601); % Time axis

preci=0.01; % If we setup the number of generated components to limit 
risl=0.05;  % the errorf

t = linspace(0,100,1601);
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

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Generation de nouvelles reponses impulsionnelles simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generation');
nsim=10;
h=zeros(nsim,length(x));

% Nombre de composantes necessaires d'aprris les mathematiques.
%h=exp(-0.5*t.^2);
N0=floor((length(x)/2)^(1/alpha)/preci/risl)+1;

vt=zeros(nsim,N0);
amp0t=zeros(nsim,N0);

% a constant depending only on alpha
w=zeros(1,length(x)-1);
for i= 1:(length(x)-1)
    w(i)=(x(i+1)-x(i))*max([mu(i),mu(i+1)]);
end

mE=sum(w);
C=(mE*(1-alpha)/(2^(alpha/2)*gamma(1+alpha/2)*gamma(2-alpha)*cos(pi*alpha/2)))^(1/alpha);

%  Generation of delays and complex amplitude.
for k=1:nsim
    disp([' --- ' int2str(nsim) ' ---  Generation - ' int2str(k)])
    vt(k,:)=Tcom_O9O9_Gener_IR2(x,N0,w);       % Delays
    amp0t(k,:)=Tcom_O9O9_Gener_RI3(N0,alpha);  % Amplitude
end

%   Generation of IR simulation
amp_sim=zeros(nsim,1601);
for k=1:nsim
    v=vt(k,:);
    amp0=amp0t(k,:);
    hg=zeros(size(x));
    
    % The complete generated impulse response
    dt = x(2)-x(1);
    for i=1:length(x)-1,
        u=find(v>(x(i)-dt)& v<(x(i)+dt));
        AA=amp0(u).*sinc((v(u)-x(i))*pi/2*dt);     %%%%  window
        hg(i)=sum(AA);
    end
    amp_sim(k,:)=hg;
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       charger les reponses impulsionnelle  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% load fichier

chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';

amp_mesure=zeros(length(pos),1601);
for j = pos
    fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(j) 'freq.mat'];
    eval(['load ' fichier]);
    H=10.^(H_amp/20).*exp(sqrt(-1)*H_phase*pi/180);
    chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
    [h,tau]=Freq2Time(H,f,2,8);

    figure
    plot(tau,abs(h))
    hold on
    title('h apres TF');
    xlabel('time ns');
    hold off
    
    amp_mesure(j-pos(1)+1,:)=h(1:1601);
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       normalisation  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amp_sim_max=zeros(1,nsim);
for i=1:nsim
    amp_sim_max=max(abs(amp0t(i,:)));
end
amp_sim_max=mean(amp_sim_max);

amp_mesure_max=zeros(1,length(pos));

for i=1:length(pos)
    amp_mesure_max=max(abs(amp_mesure(i,:)));
end
amp_mesure_max=mean(amp_mesure_max);

amp0t=amp0t/(amp_sim_max/amp_mesure_max);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       comparaison  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_set=5:5:30000;
covar=zeros(1,length(plot_set));

indice =1;
for n=plot_set
    disp([int2str(n) '/' int2str(plot_set(end))]);
    for k=1:nsim
        v=vt(k,1:n);
        amp0=amp0t(k,1:n);
        hg=zeros(size(x));
        % The complete generated impulse response
        dt = x(2)-x(1);
        for i=1:length(x)-1,
            u=find(v>(x(i)-dt)& v<(x(i)+dt));
            AA=amp0(u).*sinc((v(u)-x(i))*pi/2*dt);     %%%%  windows
            hg(i)=sum(AA);
        end
        
        for j=1:length(pos)
            covar(indice) = covar(indice) + corr2(abs(amp_mesure(j,:)),abs(hg));
        end
        covar(indice)=covar(indice)/(nsim*length(pos));
    end
    indice=indice+1;
end