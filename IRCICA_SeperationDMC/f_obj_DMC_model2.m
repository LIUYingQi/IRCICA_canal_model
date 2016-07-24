%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ce function est dans le but de calculer object function 
%%%  pour gradient descendre pour chercher les parametre optimise
%%%   dans le model DMC
%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = f_obj_DMC_model2(x,curve_h_dB,LOS)
dmc=DMC_model_v2(x(1),x(2),x(3),x(4),LOS,x(5));
obj = 0;
len = floor(0.95*length(dmc));
dmc=dmc(1,1:len);
obj = rms(curve_h_dB(1:len) -dmc );

