0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fichier est dans le but de visualiser simulation
%%% et mesure deja fait.
%%% on va comparer simulation avec N0 pas grand et mesure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  load fichier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
load ('simulation_IR_result.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  generation des reponses impulsionelles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('generation de reponse impulsionelle en changant N0 ')
dense_fh=mean(fh,1);
omega=[1:0.0125:3];       
Hws=[];
f=1:0.00125:3;
N0=1:1:5000;
amp=zeros(length(N0),length(x));
AMP=zeros(length(N0),length(f));

w=zeros(1,length(x)-1);
for i= 1:(length(x)-1)
    w(i)=(x(i+1)-x(i))*max([dense_fh(i),dense_fh(i+1)]);
end
mE=sum(w);

%%% C_alpha constant en formule 26
C=(mE*(1-alpha)/(2^(alpha/2)*gamma(1+alpha/2)*gamma(2-alpha)*cos(pi*alpha/2)))^(1/alpha);

%%% generer implusion reponse
v=generer_v(x,N0(length(N0)),w);       % Instants
h_const=generer_h_const(N0(length(N0)),alpha);   % Amplitude temp

for j=N0
    disp(['Ok - ' int2str(round(j/(N0(2)-N0(1)))+1) ' / ' int2str(length(N0))])
    %%% generer implusion reponse
    v_temp=v(1:j);       % Instants
    h_const_temp=h_const(1:j);   % amp instant
    amp(round(j/(N0(2)-N0(1)))+1,:)=generer_h(x,v_temp,h_const_temp,C);
    u=(ones(length(f),1)*amp(round(j/(N0(2)-N0(1)))+1,:)).*exp(sqrt(-1)*f'*x);
    AMP(round(j/(N0(2)-N0(1)))+1,:)=C*sum(u.');
    AMP(round(j/(N0(2)-N0(1)))+1,:)=AMP(round(j/(N0(2)-N0(1)))+1,:)/sqrt(sum(abs(AMP(round(j/(N0(2)-N0(1)))+1,:)).^2));    
    amp(round(j/(N0(2)-N0(1)))+1,:)=abs(amp(round(j/(N0(2)-N0(1)))+1,:));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = N0
    hold on  
    title(['fig - N0: ' int2str(k) ' / ' int2str(N0(length(N0)))]);
    xlabel('ns');
    ylabel('h');
	plot((abs(amp(round(k/(N0(2)-N0(1)))+1,:))))
	M(k) = getframe;
    hold off
end

figure
title('N0');
xlabel('ns');
ylabel('h');
movie(M,5)
