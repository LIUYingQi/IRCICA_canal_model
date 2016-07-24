function impuls=generer_h(pas,v,amp,C)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cette fonction est pour generer simulation pour reponse impulsionelle
%%% formule 28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=pas(2)-pas(1);
u=zeros(size(pas));
u(1:2)=0;
for i=1:length(pas)
    t0=pas(i);
    t1=t0-dt;
    t2=t0+dt;

    indice=find(v>=t1 & v<t2);
    if isempty(indice)==0
        AA=amp(indice).*sinc((v(indice)-t0)*pi/(2*dt));
        u(i)=sum(AA);
    end
end
impuls=2*pi*C*u;