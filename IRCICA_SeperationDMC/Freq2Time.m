function [h,tau]=Freq2Time(Hc,f,fenetre,surech);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ce programme permet de calculer la reponse impulsionelle du canal en 
% partant de la fonction de transfert mesure.
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

if nargin==3
   surech=1;
elseif nargin==2
   surech=1;
   fenetre=1;
end

switch fenetre
    case 1
        fen=hanning(1601)';
    case 2
        fen=boxcar(1601)';
    case 3
        fen=hamming(1601)';
    case 4
        fen=barthannwin(1601)';
    case 5
        fen=bartlett(1601)';
    case 6
        fen=blackman(1601)';
    case 7
        fen=blackmanharris(1601)';
    case 8
        fen=chebwin(1601)';
    case 9
        fen=kaiser(1601)';
end

L1=length(Hc);
Mil=floor(L1/2);
L2=L1*surech;
pas=f(2)-f(1);
pasT=1/(f(end)-f(1))/surech;
Tmax=1/pas+pasT*(surech-1);
tau=0:pasT:Tmax;
Hf=Hc.*fen;
Hf_r=[Hf(Mil+1:L1) zeros(1,(surech-1)*L1) Hf(1:Mil)];
h=(1/(tau(2)-tau(1)))*ifft(Hf_r);


