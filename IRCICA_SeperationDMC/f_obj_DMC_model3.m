%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ce function est dans le but de calculer object function 
%%%  pour gradient descendre pour chercher les parametre optimise
%%%   dans le model DMC
%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = f_obj_DMC_model3(x,curve_h_dB,LOS)
dmc=DMC_model_v3(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),LOS,x(11));
obj = 0;
len = floor(0.95*length(dmc));
dmc=dmc(1,1:len);
obj = rms(curve_h_dB(1:len) -dmc );
