%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ce function est dans le but de calculer object function 
%%%  pour gradient descendre pour chercher les parametre optimise
%%%   dans le model DMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = f_obj_DMC_chaque_mesure_v3_result(x,curve_h_dB,i)
dmc=DMC_model_v3(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12));
obj = 0;
for j=1:250
    obj = obj + rms(curve_h_dB((i-1)*250+j) -dmc );
end