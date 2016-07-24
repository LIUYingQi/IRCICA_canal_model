function f=Tcom_O9O9_2_ctes(k,N)
f=ones(1,N);
for i=2:k
    f=[f;cumsum(f(i-1,:))];
end