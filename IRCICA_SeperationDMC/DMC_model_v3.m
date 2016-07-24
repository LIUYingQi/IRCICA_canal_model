function dmc=DMC_model_v3(a0_dB,a1_dB,a2_dB,a3_dB,a4_dB,a5_dB,alpha_d,beta_d,gamma_d,theta_d,tau1_d,tau2_d)
%%%  mathmatic model version2 pour MDC
x=1:1601;
dmc = zeros(size(x)); % Preallocating enough memory for y
for i = 1:length(x)
    if  i <= length(x)*tau1_d
        dmc(i) = a0_dB+(a1_dB-a0_dB)*exp(alpha_d*(x(i)-x(length(x))*tau1_d))+a3_dB+(a4_dB-a3_dB)*exp(gamma_d*(x(i)-x(length(x))*tau2_d));
    elseif i> length(x)*tau1_d & i<= length(x)*tau2_d
        dmc(i) = a2_dB+(a1_dB-a2_dB)*exp(-beta_d*(x(i)-x(length(x))*tau1_d))+a3_dB+(a4_dB-a3_dB)*exp(gamma_d*(x(i)-x(length(x))*tau2_d));
    else
        dmc(i) = a2_dB+(a1_dB-a2_dB)*exp(-beta_d*(x(i)-x(length(x))*tau1_d))+a5_dB+(a4_dB-a5_dB)*exp(-theta_d*(x(i)-x(length(x))*tau2_d));
    end
end
