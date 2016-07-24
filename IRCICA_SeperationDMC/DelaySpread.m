function [DS2,DSp]=DelaySpread(h,t,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ce fonction est dans le but de generer delay spread de h
%%% delay spread est un propriete pour un canal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=sum(t.^p.*abs(h).^p);
b=sum(t.*abs(h))/length(t);
c=sum(abs(h).^p);
DSp=((a-b)/c)^(1/p);


p0=2;
a=sum(t.^p0.*abs(h).^p0);
b=sum(t.*abs(h).^p0);
c=sum(abs(h).^p0);
DS2=(a/c-(b/c)^p0)^(1/p0);

