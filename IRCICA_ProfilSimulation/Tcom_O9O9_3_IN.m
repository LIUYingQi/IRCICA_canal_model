function s=Tcom_O9O9_3_IN(X,lambda,h,N,alpha,m,tau)

e=exp(sqrt(-1)*(((tau*m))'*(2*pi*lambda)));
s=(abs((tau^(1/alpha))*real(((h.*X)*e))));