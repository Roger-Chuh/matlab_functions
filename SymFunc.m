function SymFunc()
syms R_c0cn R_c0cm R_cmpj R_p0pj R_p0pi t_c0cn t_c0cm t_cmpj t_p0pj t_p0pi R_cnpi t_cnpi T_cn_pi


 T_cnpi = [R_c0cn' -R_c0cn' * t_c0cn;0 1] * [R_c0cm t_c0cm;0 1] * [R_cmpj t_cmpj;0 1] * [R_p0pj' -R_p0pj'*t_p0pj;0 1] * [R_p0pi t_p0pi;0 1];

 
 R_cnpi = T_cnpi(1,1);
 t_cnpi = T_cnpi(1,2);
 
 
end