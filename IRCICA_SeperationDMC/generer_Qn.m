function f=generer_Qn(k,N)
f=ones(1,N);
for i=2:k
    f=[f;cumsum(f(i-1,:))];
end