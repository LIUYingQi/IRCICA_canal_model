function [f,H]=loadH(pos,emplacement)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   ce fonction est dans le but de load les reponse frequentielle de 
%%%   reponse impulsionnelle  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance=[2.5,2.5,2.5,3,3,3,3.5,3.5,3.5,4,4,4,4,4.5,4.5,4.5,5,5,5,5.5,5.5,6,6,6.5,6.5,7]-2.5;
% distance de mouvement pendant mesure (emplacement)
chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(pos) 'freq.mat'];
if exist(fichier)==0
    f=[];
    H=[];
    disp(['le fichier n''existe pas'])
else
    eval(['load ' fichier]);
    if (length(H_amp)==1601 && length(f)==1601)
        H=10.^(H_amp/20).*exp(sqrt(-1)*H_phase*pi/180).*exp(sqrt(-1)*(f)*2*pi*(distance(emplacement))/0.3/cos(20*pi/180));
    else
        f=[];
        H=[];
        disp(['Problem de taille'])
    end
end
