function h=Tcom_O9O9_Gener_RI3(n,alp)

%rand('state',sum(100*clock))
tmp=-log(rand(1,n));
b=(cumsum(tmp)).^(-1/alp);% the arrival times of poisson process  
                          % distributed as gamma of parameter i
h=b.*(randn(1,n)+sqrt(-1)*randn(1,n));%t2';%uu;  %iabl on the unit circle
