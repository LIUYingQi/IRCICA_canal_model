function [f,H]=Tcom_O9O9_loadH(pos,emplacement)

distance=[2.5,2.5,2.5,3,3,3,3.5,3.5,3.5,4,4,4,4,4.5,4.5,4.5,5,5,5,5.5,5.5,6,6,6.5,6.5,7]-2.5;
%chemin='E:\Donnees\P3_4_16_Circ_Droite';
%chemin='D:\';
%chemin='C:\Donnees\salle324\Lin_1_4';
% chemin='D:\Donnees\Lin_1_4';
chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';

%chemin='C:\Lin_1_4';
fichier=[chemin '\Emplacement' int2str(emplacement) '\C' int2str(pos) 'freq.mat'];
if exist(fichier,'file')==0
    f=[];
    H=[];
    disp('le fichier n''existe pas')
else
    eval(['load ' fichier]);
    if (length(H_amp)==1601 && length(f)==1601)
        H=10.^(H_amp/20).*exp(sqrt(-1)*H_phase*pi/180).*exp(sqrt(-1)*(f)*2*pi*(distance(emplacement))/0.3/cos(20*pi/180));
    else
        f=[];
        H=[];
        disp('Pb de taille')
    end
end
