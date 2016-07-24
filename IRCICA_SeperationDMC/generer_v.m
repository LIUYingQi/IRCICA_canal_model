function f=generer_v(x,N0,w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generation une simulation de reponse impulsionnel
%%% parametre v en formule 25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=w/sum(w);
I=zeros(1,N0);
for i=1:N0
    a=1;
    u=rand;
    P=w(1);
    while(u>P && a<length(w))
        a=a+1;
        P=P+w(a);
    end
    I(i)=a;
end    
f=x(I)+(x(I+1)-x(I)).*rand(1,N0);
