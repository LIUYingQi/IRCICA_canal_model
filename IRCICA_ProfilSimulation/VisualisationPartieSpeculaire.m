%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           ce ficher est dans le but de generer des simulation
%%%           pour voir l'influence de threshold N0 pour partie speculaire
%%%           et partie dense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      les formules indique dans ce programmme sont dans le papier
%%%      <<Statistical Channel Model Based on alpha-Stable Random Process
%%%      and Application tothe 60 GHz Ultra Wide Band Channel>>
%%%      Nourddine Azzaoui et laurent Clavier, Menber IEEE, IEEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

napp=30;     % Number of impuls responses used for estimating mu.
Ntest=250;  % Number of available impulse responses 

alpha=1.83; % Chosen alpha - should be estimated...
p=alpha/2.4;
x=linspace(0,100,1601); % Time axis

preci=0.01; % If we setup the number of generated components to limit 
risl=0.05;  % the errorf

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Estimation of mu.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('IndiceFH.mat');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Generation de nouvelles reponses impulsionnelles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generation');
nsim=1;
h=zeros(nsim,length(x));

% Nombre de composantes necessaires d'aprris les mathematiques.
%h=exp(-0.5*t.^2);
N0=floor((length(x)/2)^(1/alpha)/preci/risl)+1;

vt=zeros(nsim,N0);
amp0t=zeros(nsim,N0);

% a constant depending only on alpha
w=zeros(1,length(x)-1);
for i= 1:(length(x)-1)
    w(i)=(x(i+1)-x(i))*max([mu(i),mu(i+1)]);
end
mE=sum(w);
C=(mE*(1-alpha)/(2^(alpha/2)*gamma(1+alpha/2)*gamma(2-alpha)*cos(pi*alpha/2)))^(1/alpha);

%  Generation of delays and complex amplitude.
for k=1:nsim
    disp([' --- ' int2str(nsim) ' ---  Generation - ' int2str(k)])
    vt(k,:)=Tcom_O9O9_Gener_IR2(x,N0,w);       % Delays
    amp0t(k,:)=Tcom_O9O9_Gener_RI3(N0,alpha);  % Amplitude
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Test sur la premiere reponse generer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=vt(1,:);
amp0=amp0t(1,:);

hg=zeros(size(x));

% The complete generated impulse response
dt = x(2)-x(1);
for i=1:length(x)-1,
    u=find(v>(x(i)-dt)& v<(x(i)+dt));
    AA=amp0(u).*sinc((v(u)-x(i))*pi/2*dt);     %%%%  windows
    hg(i)=sum(AA);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For a different number of elements taling into account we plot the
% specular part and the dense part as well as the total RI.

dtest=[10 50 100 150 250 500];

for id=1:length(dtest);
    d0=dtest(id);
    vS=v(1:d0);     % For the specular part;
    vD=v(d0+1:N0);  % For the dense part;
    hhS=zeros(size(x));
    hhD=zeros(size(x));
    
    for i=1:length(x)-1
        uS=find(vS>x(i)-dt & vS<x(i)+dt);
        AA=amp0(uS).*sinc((v(uS)-x(i))*pi/2*dt); 
        hhS(i)=sum(AA);
        uD=find(vD>(x(i)-dt)& vD<(x(i)+dt));
        AA=amp0(d0+uD).*sinc((v(d0+uD)-x(i))*pi/2*dt);
        hhD(i)=sum(AA);
    end
    
    %     X1=20*log10(abs(hg));
    %     X2=20*log10(abs(hhS));
    %     X3=20*log10(abs(hhD));
    X1=abs(hg);
    X2=abs(hhS);
    X3=abs(hhD);
    
    figure
    subplot(311)
    plot(x,X1)
    My=max(X1);
    
    my=min(X1);
    axis([x(1) x(end) my My])
    subplot(312)
    plot(x,X2)
    title(['Specular,' int2str(d0) ' components'])
    axis([x(1) x(end) my My])
    subplot(313)
    plot(x,X3,'r')
    title('Dense')
    axis([x(1) x(end) my My]) 
end
