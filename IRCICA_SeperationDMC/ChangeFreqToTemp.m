%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ce programme permet de calculer la reponse impulsionelle du canal en 
% partant de la fonction de transfert mesure. 
% soit de domaine frequntielle a domaine temporelle
% Parametres:
%     fenetre : 0 pour une fenetre carre,
%               1 pour une fenetre de Hanning.
%     surech  : facteur de sur-ehantillonnage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hc domaine frequentielle 
% f frequence 1:1601:3   represente 57:59 GHz
% fenetre hanning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

emplacement = 14;
pos = 250;
fen = boxcar(1601);
%%% load fichier
chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(pos) 'freq.mat'];
if exist(fichier)==0
    disp(['le fichier n''existe pas'])
else
    eval(['load ' fichier]);
    if (length(H_amp)==1601 && length(f)==1601)
        disp(['load H ok'])
        H=10.^(H_amp/20).*exp(sqrt(-1)*H_phase*pi/180);
    else
        disp(['Problem de taille H'])
    end
end

chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
fichier=[chemin '\emplacement' int2str(emplacement) '\C' int2str(pos) 'temp.mat'];
if exist(fichier)==0
    disp(['le fichier n''existe pas'])
else
    eval(['load ' fichier]);
    if (length(h_amp)==1601 && length(t)==1601)
        disp(['load h ok'])
    else
        disp(['Problem de taille h'])
    end
end

%%% changement freq a temp
[h,tau]=Freq2Time(H,f,1,8);
h_amp=10.^(h_amp/20);

figure 
plot(t,h_amp,tau,abs(h))

figure 
plot(t,h_amp,t,abs(h(1:1601)))

%%% changement temp a freq

L1=length(H);
Mil=floor(L1/2);
L2=L1*8;
pas=f(2)-f(1);
pasT=1/(f(end)-f(1))/8;
Tmax=1/pas+pasT*(8-1);
tau=0:pasT:Tmax;
h = (tau(2)-tau(1))*h;
Hf = fft(h);
Hf = fftshift(Hf);
Hf = Hf(floor(1601*3.5)+1:floor(1601*4.5));
e = 0.0005;
Hf = abs(Hf).*((ones(1601,1)+e)./(fen+e))';

%%% deux figure pour voir le resultat de TF TF inverse

figure
plot(abs(Hf))

figure 
plot(10.^(H_amp/20))