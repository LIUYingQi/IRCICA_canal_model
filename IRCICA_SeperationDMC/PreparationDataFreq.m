%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de transformer FFT de freq a temp
%%% pour les 5600 mesure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% load H p to h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
H=zeros(6500,1601);
h=zeros(6500,1601*8);
f=linspace(1,3,1601);

for emplacement=1:26
    for pos=1:250
        fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(pos) 'freq.mat'];
        eval(['load ' fichier]);
        disp([ 'OK : ' int2str((emplacement-1)*250+pos) '/5600']);
        H((emplacement-1)*250+pos,:)=10.^(H_amp/20).*exp(sqrt(-1)*H_phase*pi/180);
        [h((emplacement-1)*250+pos,:),tau]=Freq2Time(H((emplacement-1)*250+pos,:),f,1,8);
    end
end

save response_temp_FFT.mat h