%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ce function est dans le but de calculer object function 
%%%  pour gradient descendre pour chercher les parametre optimise
%%%   dans le model DMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = f_obj_DMC_Richter(x,curve,LOS)
dmc=DMC_model_v1(x(1),x(2),x(3),LOS);
obj = 0;
len = floor(0.95*length(dmc));
dmc=dmc(1,1:len);
obj = rms(curve(1:len) -dmc );