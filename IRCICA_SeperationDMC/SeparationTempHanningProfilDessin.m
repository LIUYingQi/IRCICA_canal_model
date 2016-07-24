%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de voir le profil 
%%% puissance moyen et delay spread pour chaque emplacement
%%% dans la fin de ce programme on va voir le profil 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      <<Statistical Channel Model Based on alpha-Stable Random Process
%%%      and Application tothe 60 GHz Ultra Wide Band Channel>>
%%%      Nourddine Azzaoui et Laurent Clavier, Menber IEEE, IEEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% change ici pour voir i eme emplacement parmi les 26 emplacements dans la these
i=10;
%%% change ici pour definir seuil SNR
dB_seuil= 3;
%%% change ici pour definir seuil bruit
dB_bruit= -65;

%%%   load fichier result dans traitement_h.m 
tsim=linspace(0,100,1601); 
load('curve_h_dB');

DMC3 = zeros(6500,12);
h_mesure = zeros(250,1601);
h_mesure_dB = zeros(250,1601);
h_mesure_phase = zeros(250,1601);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  chargement donnes mesure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chemin='C:\Users\liuyingqi\Documents\MATLAB\E-4p_R-16';
for j=1:250
    fichier=[chemin '\Emplacement' int2str(i) '\C' int2str(j) 'temp.mat'];
    if exist(fichier)==0
        disp(['le fichier n''existe pas'])
    else
        eval(['load ' fichier]);
        if (length(h_amp)==1601 && length(t)==1601)
            disp(['chargememnt OK'])
            h_mesure_dB(j,:)=h_amp;
            h_amp=10.^(h_amp/20);
            h_mesure(j,:)=h_amp.*exp(sqrt(-1)*h_phase);
            h_mesure_phase(j,:)=h_phase;
        else
            disp(['Problem de taille'])
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   chercher DMC pour les 250 mesures  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour chaque 250 mesure
for j=1:250
    DMC3((i-1)*250+j,:)=fminsearch(@(x) f_obj_DMC_chaque_mesure_v3(x,curve_h_dB,(i-1)*250+j),[-70,-20,-40,-30,-20,-50,0.08,0.01,0.05,0.01,0.2,0.7],optimset('MaxFunEvals',600000,'MaxIter',600000));
end
% moyen de ce 250 mesure pour trouver un DMC parametrique
DMC3_result=fminsearch(@(x) f_obj_DMC_chaque_mesure_v3(x,curve_h_dB,(i-1)*250+j),[-70,-20,-40,-30,-20,-50,0.08,0.01,0.05,0.01,0.2,0.7],optimset('MaxFunEvals',600000,'MaxIter',600000));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% processus de traitement pour chercher des partie speculaire
%%% pour les 250 mesure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

power_spec_table = zeros(1,250);
power_DMC_table = zeros(1,250);
DS_table = zeros(1,250);
alpha=1.83; %1.83 alpha 
p=alpha/2.4;

for j=1:250
    DMC_h = DMC_model_v3(DMC3_result(1),DMC3_result(2),DMC3_result(3),DMC3_result(4),DMC3_result(5),DMC3_result(6),DMC3_result(7),DMC3_result(8),DMC3_result(9),DMC3_result(10),DMC3_result(11),DMC3_result(12));
    data = detrend(h_mesure_dB(j,:)-DMC_h);
    [pks,locs] = findpeaks(data,tsim,'MinPeakHeight',dB_seuil,'MinPeakDistance',0.5);
    
    %%%   chercher partie principale & enlever partie DMC
    [pks_h,locs_h] = findpeaks(h_mesure_dB(j,:),tsim,'MinPeakDistance',0.5);

    %%%   enlever partie avant pic principal
    max_pks = pks_h(1);
    for k=1:length(pks_h)
        if pks_h(k)>max_pks
            max_pks = pks_h(k);
            max_loc = locs_h(k);
        end
    end
    k=find(abs(locs-max_loc) < 0.1);
    pks = pks(k:length(pks));
    locs = locs(k:length(locs));

    %%%   enlever partie bruit
    len=length(locs);
    delete_h=[];
    for k=1:len
        if h_mesure_dB(j,round(1601*locs(k)/100)) < dB_bruit
            delete_h=[delete_h k];
        end
    end
    locs(delete_h)=[];
    pks(delete_h)=[];

    %%% profil puissance et delay spread
    a_t=h_mesure(j,round(locs/100*1601))-DMC_h(round(locs/100*1601));
    tau_t=locs;
    theta_t=h_phase(round(locs/100*1601));
    G_w=hanning(length(t));
    G_w=[G_w' zeros(1,7*length(t))];
    g_t=ifft(G_w);
    g_t=ifftshift(g_t);
    g_t=(length(t)/100)*g_t;
    partie_spectrale_real=zeros(1,length(t));
    partie_spectrale_img=zeros(1,length(t));
    partie_spectrale = complex(partie_spectrale_real,partie_spectrale_img);
    for k=1:length(locs)
        partie_spectrale(round(tau_t(k)/100*1601))=a_t(k).*exp(sqrt(-1)*theta_t(k)*pi/180);
    end
    partie_spectrale=conv(partie_spectrale,g_t,'same');
    DMC_h=h_mesure(j,:)-abs(partie_spectrale);
    
    %%% resultat
    power_spec = mean(abs(partie_spectrale).^2);
    power_DMC = mean(DMC_h.^2);
    DS = DelaySpread(h_mesure(j,:),tsim,p);
    
    power_spec_table(j) = power_spec;
    power_DMC_table(j) = power_DMC;
    DS_table(j) = DS;
end

power_spec = mean(power_spec_table);
power_DMC = mean(power_DMC_table);
DS = mean(DS_table);
