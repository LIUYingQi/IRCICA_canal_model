%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de separer signal avec boxcar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% change ici pour voir i eme mesure   
num_mesure=3500;
emplacement = floor(num_mesure/250)+1;    
pos = mod(num_mesure,250)+1;    
%%% change ici pour definir seuil SNR
dB_seuil= 4;
%%% change ici pour definir seuil bruit
dB_bruit= -70;
%%% definir fenetre utiliser (1,hanning 2,boxcar)
fen = 2;
fens = boxcar(1601);
e = 0.0005;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  restaurer signal initial sans influence de fenetre hanning
%%%  avec 2 partie sepsrer par facon 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% load fichier

chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(pos) 'freq.mat'];
eval(['load ' fichier]);
H=10.^(H_amp/20).*exp(sqrt(-1)*H_phase*pi/180);
chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(pos) 'temp.mat'];
eval(['load ' fichier]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% changement domaine freq a domaine temp (pas de hanning mais boxcar)
%%% h = h_speculaire + h_dense
%%% calculer h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,tau]=Freq2Time(H,f,fen,8);
h_amp=10.^(h_amp/20);
figure 
plot(t,h_amp,tau,abs(h))
hold on
title('h apres TF');
xlabel('time ns');
legend('h-amp dans temp.mat','h apres inverse TF');
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculer h_speculaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_amp = abs(h(1:length(t)));
h_phase = angle(h(1:length(t)));

%%% model DMC et curve-fitting

h=h_amp.*exp(sqrt(-1)*h_phase*pi/180); 
h_dB=db(h_amp);   %% load
h_amp_dB = h_dB;
remove_pike_h_dB=medfilt1(h_dB,10);   %% enlever pics attenues
smt_h_dB = sgolayfilt(remove_pike_h_dB,3,41);    %%  smooth
poly=polyfit(t,smt_h_dB,100);      %%   polyfit
curve_h_dB=polyval(poly,t);    %%  exprimer en polynome

figure
hold on
plot(t,db(abs(h)),t,curve_h_dB)
title('curve-fitting pour DMC model');
xlabel('time ns');
legend('h','curve DMC');
hold off

DMC=fminsearch(@(x) f_obj_DMC_res(x,curve_h_dB),[-40,-65,-80,0.05,0.2,0.01,],optimset('MaxFunEvals',600000,'MaxIter',600000));

figure
hold on
plot(t,curve_h_dB)
plot(t,DMC_model_v2(DMC(1),DMC(2),DMC(3),DMC(4),DMC(5),DMC(6)));
title('DMC model');
xlabel('time ns');
legend('curve','DMC');
hold off

figure
data = detrend(h_amp_dB-DMC_model_v2(DMC(1),DMC(2),DMC(3),DMC(4),DMC(5),DMC(6)));
[pks,locs,widths,proms] = findpeaks(data,t,'MinPeakHeight',dB_seuil,'MinPeakDistance',0.5);
findpeaks(data,t,'MinPeakHeight',dB_seuil,'MinPeakDistance',0.5,'Annotate','extents');
xlabel('t');
ylabel('partie speculairee');

%%%   chercher partie principale & enlever partie DMC
[pks_h,locs_h] = findpeaks(h_amp_dB,t,'MinPeakDistance',0.5);

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
plot(t,h_amp_dB)
h_amp_dB_locs=h_amp_dB(round(1601*locs_h(:)/100));
stem(locs_h,h_amp_dB_locs);
hold off

partie_DMC = DMC_model_v2(DMC(1),DMC(2),DMC(3),DMC(4),DMC(5),DMC(6));
pks_h = 10.^(pks_h/20);
partie_DMC = 10.^(partie_DMC/20);

a_t=pks_h-partie_DMC(round(locs_h/100*1601));
tau_t=locs_h;
theta_t=h_phase(round(locs_h/100*1601));
partie_spectrale_real=zeros(1,length(t));
partie_spectrale_img=zeros(1,length(t));
partie_spectrale = complex(partie_spectrale_real,partie_spectrale_img);
for j=1:length(locs_h)
    partie_spectrale(round(tau_t(j)/100*1601))=a_t(j).*exp(sqrt(-1)*theta_t(j));
end

figure
hold on
plot(t,h_amp,t,abs(partie_spectrale),t,partie_DMC)
title('partie speculaire');
xlabel('time ns');
legend('h','h_s','DMC model');
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculer h_dense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partie_spectrale = [partie_spectrale complex(zeros(1,length(t)*7),zeros(1,length(t)*7))];
L1=length(H);
Mil=floor(L1/2);
L2=L1*8;
pas=f(2)-f(1);
pasT=1/(f(end)-f(1))/8;
Tmax=1/pas+pasT*(8-1);  
tau=0:pasT:Tmax;
H_s = fft(partie_spectrale);
H_s = fftshift(H_s);
H_s = H_s(floor(1601*3.5)+1:floor(1601*4.5));
H_s = H_s.*((ones(1601,1)+e)./(fens+e))';

H_d = H - H_s;

figure
hold on
plot(f+56,abs(H))
plot(f+56,abs(H_s))
plot(f+56,abs(H_d))
title('domaine frequentielle');
xlabel('f GHz ');
legend('|H|','|H_s|','|H_d|');
hold off

[h_d,tau]=Freq2Time(H_d,f,fen,8);
partie_spectrale = partie_spectrale(1:length(t));
partie_dense = h_d(1:length(t));
h = h(1:length(t));

figure
hold on
plot(t,abs(partie_spectrale))
plot(t,abs(partie_dense))
plot(t,abs(h))
title('separation de partie speculaire et partie dense enlever fenetre hanning');
xlabel('time ns');
legend('h_s','h_d','h');
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% tableau pour voir puissance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('puissance ');
power_spec = mean(abs(partie_spectrale).^2);
power_DMC = mean(abs(partie_dense).^2);
T = table([power_spec],[power_DMC],'VariableNames',{'power_spec' 'power_DMC'})
rapport = power_spec/power_DMC

