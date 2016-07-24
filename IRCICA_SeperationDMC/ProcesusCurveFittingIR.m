%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       ce ficher est pour dessiner reponse inpulsionnelle de tous les
%%%       resultat obtenu pour voir globalement les proprietes de canal
%%%       on va dessiner plot original, apres remove pikes, apres curve fitting   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   load, calculation
for i=0:1:6499
    %load reponse impulsionnelle frequentielle
    [f,H(i+1,:)]=loadH(mod(i,250)+1,floor(i/250)+1);
    
    %load reponse impulsionelle temporelle
    chemin='C:\Users\liuyingqi\Desktop\E-4p_R-16';
    fichier=[chemin '\emplacement' int2str(floor(i/250)+1) '\C' int2str(mod(i,250)+1) 'temp.mat'];
    if exist(fichier)==0
        t=[];
        h=[];
        h_dB=[];
        Sdisp(['le fichier n''exsite pas - Emp ' int2str(floor(i/250)+1) ' Pos ' int2str(mod(i,250)+1)])
    else
        disp(['OK - Emp ' int2str(floor(i/250)+1) ' Pos ' int2str(mod(i,250)+1)])
        eval(['load ' fichier]);
        if (length(h_amp)==1601 && length(f)==1601)
            h(i+1,:)=10.^(h_amp/20).*exp(sqrt(-1)*h_phase*pi/180); 
            h_dB(i+1,:)=h_amp;   %% load
            remove_pike_h_dB(i+1,:)=medfilt1(h_dB(i+1,:),10);   %% enlever pics attenues
            smt_h_dB(i+1,:) = sgolayfilt(remove_pike_h_dB(i+1,:),3,41);    %%  smooth
            poly=polyfit(tsim,smt_h_dB(i+1,:),100);      %%   polyfit
            curve_h_dB(i+1,:)=polyval(poly,tsim);    %%  exprimer en polynome
        else
            disp(['Pb de taille - Emp ' int2str(floor(i/250)+1) ' Pos ' int2str(mod(i,250)+1)])
        end
    end
end

%%      plot

i=500;  % voir i eme data mesure dans les 26*250 mesure 

subplot(4,1,1)
plot(tsim,h_dB(i,:))
axis([0 100 -120 -20])
subplot(4,1,2)
plot(tsim,remove_pike_h_dB(i,:))
axis([0 100 -120 -20])
subplot(4,1,3)
plot(tsim,smt_h_dB(i,:))
axis([0 100 -120 -20])
subplot(4,1,4)
plot(tsim,curve_h_dB(i,:))
axis([0 100 -120 -20])