%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de enlever prtie DMC (partie bruit colore, partie difficile a differer)
%%% pour voir si notre DMC modele est bon est si les
%%% pics peut etre identifie apres enlever DMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% change ici pour voir i eme mesure   
num_mesure=3500;
emplacement = floor(num_mesure/250)+1;    
pos = mod(num_mesure,250)+1;    
%%% change ici pour definir seuil SNR
dB_seuil= 3;
%%% change ici pour definir seuil bruit
dB_bruit= -70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  traitement du signal pour smooth et pour curve-fitting 
%%%  les procesus sont deja fait dans taitement_h.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   load fichier result dans traitement_h.m 
load('simulation_IR_result');
load('curve_h_dB');
load('response_temp_dB');
load('response_temp');
load('DMC_parametre.mat');

%%%  chargement donnes mesure
chemin='C:\Users\liuyingqi\Documents\MATLAB\E-4p_R-16';
fichier=[chemin '\Emplacement' int2str(floor(num_mesure/250)+1) '\C' int2str(mod(num_mesure,250)+1) 'temp.mat'];
if exist(fichier)==0
    disp(['le fichier n''existe pas'])
else
    eval(['load ' fichier]);
    if (length(h_amp)==1601 && length(t)==1601)
        disp(['chargememnt OK'])
        h_amp_dB=h_amp;
        h_amp=10.^(h_amp/20);
        h=h_amp.*exp(sqrt(-1)*h_phase);
    else
        t=[];
        h_amp=[];
        h_phase=[];
        disp(['Problem de taille'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% touver un modele pour DMC
%%% on propose 3 types de DMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   optimization   %%%
DMC(num_mesure,:)=fminsearch(@(x) f_obj_DMC_chaque_mesure_v1(x,curve_h_dB,num_mesure),[-100,-40,0.003,0.2]);
DMC2(num_mesure,:)=fminsearch(@(x) f_obj_DMC_chaque_mesure_v2(x,curve_h_dB,num_mesure),[-70,-30,-90,0.005,0.2,0.05],optimset('MaxFunEvals',600000,'MaxIter',600000));
DMC3(num_mesure,:)=fminsearch(@(x) f_obj_DMC_chaque_mesure_v3(x,curve_h_dB,num_mesure),[-70,-20,-40,-30,-20,-50,0.08,0.01,0.05,0.01,0.2,0.7],optimset('MaxFunEvals',600000,'MaxIter',600000));

%%%    plot  h_amp_dB & les 3 DMC model  %%%
figure
subplot(3,2,1)
plot(tsim,h_amp_dB)
title('h_amp_dB')
xlabel('t')
ylabel('dB')
axis([0 100 -120 -20])
subplot(3,2,2)
hold on
plot(tsim,DMC_model_v1(DMC(num_mesure,1),DMC(num_mesure,2),DMC(num_mesure,3),DMC(num_mesure,4)))
plot(tsim,DMC_model_v2(DMC2(num_mesure,1),DMC2(num_mesure,2),DMC2(num_mesure,3),DMC2(num_mesure,4),DMC2(num_mesure,5),DMC2(num_mesure,6)))
plot(tsim,DMC_model_v3(DMC3(num_mesure,1),DMC3(num_mesure,2),DMC3(num_mesure,3),DMC3(num_mesure,4),DMC3(num_mesure,5),DMC3(num_mesure,6),DMC3(num_mesure,7),DMC3(num_mesure,8),DMC3(num_mesure,9),DMC3(num_mesure,10),DMC3(num_mesure,11),DMC3(num_mesure,12)))
title('3 model');
xlabel('t');
ylabel('dB');
legend('model1','model2','model3');
axis([0 100 -120 -20])
subplot(3,2,3)
plot(tsim,h_amp_dB-DMC_model_v1(DMC(num_mesure,1),DMC(num_mesure,2),DMC(num_mesure,3),DMC(num_mesure,4)),tsim,detrend(h_amp_dB-DMC_model_v1(DMC(num_mesure,1),DMC(num_mesure,2),DMC(num_mesure,3),DMC(num_mesure,4))))
axis([0 100 -30 30])
title('model1 enlever DMC')
xlabel('t')
ylabel('dB')
legend('enlever DMC','detrend');
subplot(3,2,4)
plot(tsim,h_amp_dB-DMC_model_v2(DMC2(num_mesure,1),DMC2(num_mesure,2),DMC2(num_mesure,3),DMC2(num_mesure,4),DMC2(num_mesure,5),DMC2(num_mesure,6)),tsim,detrend(h_amp_dB-DMC_model_v2(DMC2(num_mesure,1),DMC2(num_mesure,2),DMC2(num_mesure,3),DMC2(num_mesure,4),DMC2(num_mesure,5),DMC2(num_mesure,6))))
axis([0 100 -30 30])
title('model2 enlever DMC')
xlabel('t')
ylabel('dB')
legend('enlever DMC','detrend');
subplot(3,2,5)
plot(tsim,h_amp_dB-DMC_model_v3(DMC3(num_mesure,1),DMC3(num_mesure,2),DMC3(num_mesure,3),DMC3(num_mesure,4),DMC3(num_mesure,5),DMC3(num_mesure,6),DMC3(num_mesure,7),DMC3(num_mesure,8),DMC3(num_mesure,9),DMC3(num_mesure,10),DMC3(num_mesure,11),DMC3(num_mesure,12)),tsim,detrend(h_amp_dB-DMC_model_v3(DMC3(num_mesure,1),DMC3(num_mesure,2),DMC3(num_mesure,3),DMC3(num_mesure,4),DMC3(num_mesure,5),DMC3(num_mesure,6),DMC3(num_mesure,7),DMC3(num_mesure,8),DMC3(num_mesure,9),DMC3(num_mesure,10),DMC3(num_mesure,11),DMC3(num_mesure,12))))
axis([0 100 -30 30])
title('model3 enlever DMC')
xlabel('t')
ylabel('dB')
legend('enlever DMC','detrend');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% processus de traitement pour chercher des partie speculaire
%%% ici on choisit de utiliser 3 er modele DMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
data = detrend(h_amp_dB-DMC_model_v3(DMC3(num_mesure,1),DMC3(num_mesure,2),DMC3(num_mesure,3),DMC3(num_mesure,4),DMC3(num_mesure,5),DMC3(num_mesure,6),DMC3(num_mesure,7),DMC3(num_mesure,8),DMC3(num_mesure,9),DMC3(num_mesure,10),DMC3(num_mesure,11),DMC3(num_mesure,12)));
[pks,locs,widths,proms] = findpeaks(data,tsim,'MinPeakHeight',dB_seuil,'MinPeakDistance',0.5);
findpeaks(data,tsim,'MinPeakHeight',dB_seuil,'MinPeakDistance',0.5,'Annotate','extents');
xlabel('tsim');
ylabel('partie speculairee');

%%%   chercher partie principale & enlever partie DMC
[pks_h,locs_h] = findpeaks(h_amp_dB,tsim,'MinPeakDistance',0.5);

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
    if h_amp_dB(round(1601*locs_h(k)/100)) < dB_bruit
        delete_h=[delete_h k];
    end
end
locs_h(delete_h)=[];
pks_h(delete_h)=[];

%%%  enlever partie SNR pas suffisant par rapport a DMC
len=length(locs_h);
delete_h=[];
for k=1:len
    if data(round(1601*locs_h(k)/100)) < dB_seuil
        delete_h=[delete_h k];
    end
end
locs_h(delete_h)=[];
pks_h(delete_h)=[];

figure
hold on
plot(tsim,h_amp_dB)
h_amp_dB_locs=h_amp_dB(round(1601*locs_h(:)/100));
stem(locs_h,h_amp_dB_locs);
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  analyse pour la puissance dans domaine f
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% speculaire analyse pour H total et H speculaire entre 57GHz - 59GHz
[f,H]=loadH(mod(num_mesure,250)+1,floor(num_mesure/250)+1);
%%% les parametre pour signal 
t_sampling= 100e-9/1600;    % 6.25*e-11 
f_sampling= 1/t_sampling;   % 1.6*e10 Hz  16 GHz
time=tsim;                  % sampling time
L=length(tsim);             % length

%%% 
figure 
hold on
f=f+56;
plot(f,abs(H));
title('H mesure entre 57GHz et 59GHz');
xlabel('f (GHz)');
ylabel('dB');
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  separer DMC et partie speculaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 1 er methode
%%% h soustraire DMC directement comme signal et donc laisse DMC

figure
hold on
plot(tsim,abs(h_amp))   %%%h  
h_locs=abs(h_amp(round(locs_h(:)/100*1601)));
h_spec=abs(h_amp.*(ismember(t,locs_h)));
plot(t,h_spec)    %%%   h speculairee
stem(locs_h,h_locs)
DMC_h_1=db2mag(DMC_model_v3(DMC3(num_mesure,1),DMC3(num_mesure,2),DMC3(num_mesure,3),DMC3(num_mesure,4),DMC3(num_mesure,5),DMC3(num_mesure,6),DMC3(num_mesure,7),DMC3(num_mesure,8),DMC3(num_mesure,9),DMC3(num_mesure,10),DMC3(num_mesure,11),DMC3(num_mesure,12)));
plot(tsim,DMC_h_1)
title('h mesure network analyzor travailler entre 57GHz et 59GHz');
xlabel('time ns');
hold off

figure
hold on
spec_h_1=h_amp-DMC_h_1;
for j=1:length(spec_h_1);
    if spec_h_1(j)<0;
        spec_h_1(j)=0;
    end
end
plot(tsim,DMC_h_1)
plot(tsim,spec_h_1)
title('1 er facon pour separer DMC et pics');
xlabel('time ns');
legend('DMC','composant speculaire');
hold off

disp('1 er facon pour separer DMC et pics : separer directement ');
power_spec = mean(spec_h_1.^2);
power_DMC = mean(DMC_h_1.^2);
T = table([power_spec],[power_DMC],'VariableNames',{'power_spec' 'power_DMC'})
rapport = power_spec/power_DMC

%%% 2 er methode
%%% soustraire parite speculaire dans domaine f

a_t=h_locs;
tau_t=locs_h;
theta_t=h_phase(round(locs_h/100*1601));
G_w=hanning(length(t));
G_w=[G_w' zeros(1,7*length(t))];
g_t=ifft(G_w);
g_t=ifftshift(g_t);
g_t=(length(t)/100)*g_t;
partie_spectrale_real=zeros(1,length(t));
partie_spectrale_img=zeros(1,length(t));
partie_spectrale = complex(partie_spectrale_real,partie_spectrale_img);
for j=1:length(locs_h)
    partie_spectrale(round(tau_t(j)/100*1601))=a_t(j).*exp(sqrt(-1)*theta_t(j)*pi/180);
end
partie_spectrale=conv(partie_spectrale,g_t,'same');
DMC_h_2=h_amp-abs(partie_spectrale);

figure
hold on
plot(t,abs(h_amp))   %%%h  
plot(t,abs(partie_spectrale))    %%%   h speculairee
plot(tsim,DMC_h_2)
title('2 er facon pour separer DMC et pics');
xlabel('time ns');
legend('h','composant speculaire','DMC');
hold off

disp('2 er facon pour separer DMC et pics : separer directement ');
power_spec = mean(abs(partie_spectrale).^2);
power_DMC = mean(abs(h_amp).^2)-power_spec;
T = table([power_spec],[power_DMC],'VariableNames',{'power_spec' 'power_DMC'})
rapport = power_spec/power_DMC

%%% 3 er methode
%%% soustraire parite speculaire dans domaine f sans DMC

a_t=h_locs-DMC_h_1(round(locs_h/100*1601));
tau_t=locs_h;
theta_t=h_phase(round(locs_h/100*1601));
G_w=hanning(length(t));
G_w=[G_w' zeros(1,7*length(t))];
g_t=ifft(G_w);
g_t=ifftshift(g_t);
g_t=(length(t)/100)*g_t;
partie_spectrale_real=zeros(1,length(t));
partie_spectrale_img=zeros(1,length(t));
partie_spectrale = complex(partie_spectrale_real,partie_spectrale_img);
for j=1:length(locs_h)
    partie_spectrale(round(tau_t(j)/100*1601))=a_t(j).*exp(sqrt(-1)*theta_t(j)*pi/180);
end
partie_spectrale=conv(partie_spectrale,g_t,'same');
DMC_h_3=h_amp-abs(partie_spectrale);

figure
hold on
plot(t,abs(h_amp))   %%%h  
plot(t,abs(partie_spectrale))    %%%   h speculairee
plot(tsim,DMC_h_3)
title('3 er facon pour separer DMC et pics');
xlabel('time ns');
legend('h','composant speculaire','DMC');
hold off

disp('3 er facon pour separer DMC et pics : separer directement ');
power_spec = mean(abs(partie_spectrale).^2);
power_DMC = mean(abs(h_amp).^2) - power_spec;
T = table([power_spec],[power_DMC],'VariableNames',{'power_spec' 'power_DMC'})
rapport = power_spec/power_DMC

