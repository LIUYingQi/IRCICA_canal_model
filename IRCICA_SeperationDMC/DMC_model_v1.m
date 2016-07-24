function dmc=DMC_model_v1(a0_dB,a1_dB,beta_d,tau_d)
%%%  mathmatic model pour MDC
x=1:1601;
dmc = zeros(size(x)); % Preallocating enough memory for y
for i = 1:length(x)
    if  i <= length(x)*tau_d
        dmc(i) = a0_dB ;
    else
        dmc(i) = a0_dB+(a1_dB-a0_dB)*exp(-beta_d*(x(i)-x(length(x))*tau_d));
    end
end