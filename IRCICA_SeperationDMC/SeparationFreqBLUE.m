%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de enlever prtie DMC 
%%% en utilisant methode BLUE dans page 100 dans la these de Richter
%%% c'est un methode dans domaine freq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all
clc;
warning off all

%%% change ici pour separer i eme mesure
objet = 2450;

%%% change ici pour voir empalcement
emplacement = floor(objet/250)+1;  
pos = mod(objet,250)+1;    
num_emplacememnt = 26;
%%% change ici pour definir seuil SNR
dB_seuil= 0.5;
%%% change ici pour definir seuil bruit
dB_bruit= -70;

num_pos = 250;
%%% objet wfenetre pour calculer valeur moyen
if(pos<25)
    num_mesure=num_pos*(emplacement-1)+1:1:num_pos*(emplacement-1)+50;
    num_pos_mesure=pos;
elseif(pos>=225)
    num_mesure=num_pos*(emplacement-1)+201:1:num_pos*(emplacement-1)+250;
    num_pos_mesure=pos-200;
else
    num_mesure=num_pos*(emplacement-1)+pos-25:1:num_pos*(emplacement-1)+pos+24;
    num_pos_mesure=pos-(num_pos*(emplacement-1)+pos-25);
end

%%%  chargement donnes mesure
load('response_freq.mat');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% preparation de h temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=zeros(6500,1601*8);
f=linspace(1,3,1601);

%%% inverse fft pour trouver h temp
% h = ifft(H.*repmat(hann(Mf).',6500,1),[],2);
h = ifft(H,[],2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  pre - definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mf = 1601;
bw = 2e9;
df = bw/(Mf-1);
dt = 1/bw;
delay = (0:Mf-1)*dt;
f   = linspace(1e9,3e9,1601);
t = linspace(0,100,1601);
tn = delay*1e9;

% covariance matrice
R = H(num_mesure,:)'*H(num_mesure,:)/250;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  visualisation pour les 250 test RI dans un emplacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
title('emplacement RI dB')
xlabel('t/ns')
ylabel('dB')
for i=num_mesure
    plot(tn,db(h(i,1:1601)))
end
xlim([0 150]);
hold off

h_dB_mean = 10*log10(mean(abs(h(num_mesure,1:1601)).^2,1));

figure
hold on
title('emplacement RI mena dB')
xlabel('t/ns')
ylabel('dB')
plot(tn,h_dB_mean(1:1601));
hold off


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  processus dans domaine temp pour curve fitting model parametrique 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% trouver LOS (main pic) 
LOS_loc = round(find(h_dB_mean==max(h_dB_mean))*100/1601);
LOS_loc = LOS_loc/100;
%%% curve fitting
remove_pike_h_dB=medfilt1(h_dB_mean,10);   %% enlever pics attenues
smt_h_dB = sgolayfilt(remove_pike_h_dB,3,41);    %%  smooth
poly=polyfit(tn,smt_h_dB,100);      %%   polyfit
curve_h_dB_mean=polyval(poly,tn);    %%  exprimer en polynome

figure
hold on
title('emplacement RI mena dB')
xlabel('t')
ylabel('dB')
plot(tn,h_dB_mean,tn,curve_h_dB_mean);
legend('h-dB-mean','curve-fitting');
hold off

%%% 1 er DMC model
DMC_1=fminsearch(@(x) f_obj_DMC_Richter(x,curve_h_dB_mean,LOS_loc),[-90,-50,0.01],optimset('MaxFunEvals',600000,'MaxIter',600000));
figure
hold on
title('model 1')
xlabel('t')
ylabel('dB')
plot(tn,h_dB_mean,tn,DMC_model_v1(DMC_1(1),DMC_1(2),DMC_1(3),LOS_loc));
legend('h-dB-mean','model1');
hold off

%%% 2 eme DMC model
DMC_2=fminsearch(@(x) f_obj_DMC_model2(x,curve_h_dB_mean,LOS_loc),[-70,-50,-100,0.001,0.001],optimset('MaxFunEvals',600000,'MaxIter',600000));
figure
hold on
title('model 2')
xlabel('t')
ylabel('dB')
plot(tn,h_dB_mean,tn,DMC_model_v2(DMC_2(1),DMC_2(2),DMC_2(3),DMC_2(4),LOS_loc,DMC_2(5)));
legend('h-dB-mean','model1');
hold off

%%% 3 eme DMC model
DMC_3=fminsearch(@(x) f_obj_DMC_model3(x,curve_h_dB_mean,LOS_loc),[-70,-20,-40,-30,-20,-50,0.08,0.01,0.05,0.001,0.055],optimset('MaxFunEvals',600000,'MaxIter',600000));
figure
hold on
title('model 3')
xlabel('t')
ylabel('dB')
plot(tn,h_dB_mean,tn,DMC_model_v3(DMC_3(1),DMC_3(2),DMC_3(3),DMC_3(4),DMC_3(5),DMC_3(6),DMC_3(7),DMC_3(8),DMC_3(9),DMC_3(10),LOS_loc,DMC_3(11)));
legend('h-dB-mean','model1');
hold off

%%%   chercher partie principale & enlever partie DMC
data = detrend(h_dB_mean-DMC_model_v3(DMC_3(1),DMC_3(2),DMC_3(3),DMC_3(4),DMC_3(5),DMC_3(6),DMC_3(7),DMC_3(8),DMC_3(9),DMC_3(10),LOS_loc,DMC_3(11)));
[pks_h,locs_h] = findpeaks(h_dB_mean(1:201),tn(1:201),'MinPeakDistance',1);

%%%   enlever partie avant pic principal
max_pks = pks_h(1);
for j=1:length(pks_h)
    if pks_h(j)>max_pks
        max_pks = pks_h(j);
        max_loc = locs_h(j);
    end
end
j=find(abs(locs_h-max_loc) < 0.1);
pks_h = pks_h(j:length(pks_h));
locs_h = locs_h(j:length(locs_h));

%%%   enlever partie bruit
len=length(locs_h);
delete_h=[];
for k=1:len
    if h_dB_mean(round(1601*locs_h(k)/800)+1) < dB_bruit
        delete_h=[delete_h k];
    end
end
locs_h(delete_h)=[];
pks_h(delete_h)=[];

%%%  enlever partie SNR pas suffisant par rapport a DMC
len=length(locs_h);
delete_h=[];
for k=1:len
    if data(round(1601*locs_h(k)/800)+1) < dB_seuil
        delete_h=[delete_h k];
    end
end
locs_h(delete_h)=[];
pks_h(delete_h)=[];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  processus dans domaine temp pour trouver temp d'arrive cluster tau 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = locs_h;

figure
hold on
plot(tn,h_dB_mean)
h_dB_mean_locs=h_dB_mean(round(1601*locs_h(:)/800)+1);
h_dB_mean_locs=stem(locs_h,h_dB_mean_locs,'r');
xlim([0 100])
ylim([-100 0])
set(h_dB_mean_locs,'basevalue',-200)
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  processus dans domaine freq en utilisant BLUE pour enlever partie speculaire 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BLUE
tau = tau*1e-9;
B = exp(-2*pi*1i*f'*tau)/sqrt(Mf);
As  = B'*H(num_mesure,:).';
Xs = B*As;

%%% BLUE avec covariance
% AsR = (B'/R*B)\B'/R*H(num_mesure,:).';
% XsR = B*AsR;
% Hn = H - Xs';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualisation resultat BLUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  domaine temp
figure
hold on
plot(tn,h_dB_mean)
plot(tn,20*log10(abs(ifft(H(num_mesure(1),:).'))),'k')
h_dB_mean_locs=h_dB_mean(round(1601*locs_h(:)/800)+1)
stem(locs_h,h_dB_mean_locs,'k');
plot(tn,10*log10(mean(abs(ifft(Xs)).^2,2)),'r')
h = stem(tau*1e9,10*log10(mean(abs(As).^2,2)/Mf),'r');
% h = stem(tau*1e9,10*log10(abs(As(:,1)).^2/Mf),'r');
% plot(delay*1e9,20*log10(abs(ifft(XsR))),'--r')
% h = stem(tau*1e9,20*log10(abs(AsR)/sqrt(Mf)),'-g');
set(h,'basevalue',-200)
hold off
xlim([0 100])
ylim([-100 0])

%%% domaine freq
figure
hold on
title('freq domaine')
xlabel('t')
ylabel('dB')
H=H(num_mesure,:);
plot(f,20*log10(abs(H(num_pos_mesure,:))));
Xs=Xs.';
plot(f,20*log10(abs(Xs(num_pos_mesure,:))));
hold off

