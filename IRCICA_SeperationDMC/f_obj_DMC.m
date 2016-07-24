%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ce function est dans le but de calculer object function 
%%%  pour gradient descendre pour chercher les parametre optimise
%%%   dans le model DMC
%%%   x(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = f_obj_DMC(x,curve_h_dB)
dmc=DMC_model(x(1),x(2),x(3),x(4));
obj = 0;
len = floor(0.95*length(dmc));
dmc=dmc(1,1:len);
for i=1:1:6500
    obj = obj + rms(curve_h_dB(i,1:len) -dmc );
end

