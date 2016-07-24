function f=Tcom_O9O9_Gener_IR2(x,N,w)

w=w/sum(w);
I=zeros(1,N);
for i=1:N
    a=1;
    u=rand;
    P=w(1);
    while(u>P && a<length(w))
        a=a+1;
        P=P+w(a);
    end
    I(i)=a;
end    
f=x(I)+(x(I+1)-x(I)).*rand(1,N);
