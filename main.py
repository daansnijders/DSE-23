# -*- coding: utf-8 -*-
"""
Main program of DSE-23
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""




from inputs.concept_1 import *
from inputs.constants import *
from inputs.performance_inputs import *






#CLASS2 WEIGHT AND CG LOCATION

#configuration 1 





#config1_iteration_diff=50
#config3_iteration_diff=50
#config2_iteration_diff=50
#
#
#
#config1_class2     = Class2_weight(1,N_pax[0],MTOW[0],M_carried_canard_MZF[0],min(M_MZF), n_max[0],V_dive[0],M_fuel[0], max(T_req), l_f[0],d_f_inner,d_f_outer,l_cabin[0], l_h[0], S, S_c[0], b, b_c[0], S_v[0],S_h[0],Cr_t,Cr_t_c[0],lambda_2_rad,lambda_h_2_rad[0], lambda_v_2_rad[0],lambda_c_2_rad, S_fus[0])     
#
#config1_M_structural            =config1_class2.structural_mass()
#config1_M_powerplant            =config1_class2.powerplant_mass()
#config1_M_fixedeq               =config1_class2.fixed_equipment_mass()
#config1_class2_OEW              =config1_class2.OEW(config1_M_structural,config1_M_powerplant,config1_M_fixedeq)
#cg1 = get_cg(x_le_MAC,config1_class2)   
#
#config3_class2     = Class2_weight(3,N_pax[2],MTOW[2],M_carried_canard_MZF[2],min(M_MZF), n_max[2],V_dive[2],M_fuel[2], max(T_req), l_f[2],d_f_inner,d_f_outer,l_cabin[2], l_h[2], S, S_c[2], b, b_c[2], S_v[2],S_h[2],Cr_t,Cr_t_c[2],lambda_2_rad,lambda_h_2_rad[2], lambda_v_2_rad[2],lambda_c_2_rad, S_fus[2])
#
#config3_M_structural            =config3_class2.structural_mass()
#config3_M_powerplant            =config3_class2.powerplant_mass()
#config3_M_fixedeq               =config3_class2.fixed_equipment_mass()
#config3_class2_OEW              =config3_class2.OEW(config3_M_structural,config3_M_powerplant,config3_M_fixedeq)
#cg3 = get_cg(x_le_MAC,config3_class2) 
#
#config2_class2     = Class2_weight(2,N_pax[1],MTOW[1],M_carried_canard_MZF[1],min(M_MZF), n_max[1],V_dive[1],M_fuel[1], max(T_req), l_f[1],d_f_inner,d_f_outer,l_cabin[1], l_h[1], S, S_c[1], b, b_c[1], S_v[1],S_h[1],Cr_t,Cr_t_c[1],lambda_2_rad,lambda_h_2_rad[1], lambda_v_2_rad[1],lambda_c_2_rad, S_fus[1])
#config2_M_structural            =config2_class2.structural_mass()
#config2_M_powerplant            =config2_class2.powerplant_mass()
#config2_M_fixedeq               =config3_M_fixedeq                                                 #needed as the fixed equipment needs to be the same for both
#config2_class2_OEW              =config2_class2.OEW(config2_M_structural,config2_M_powerplant,config2_M_fixedeq)
#cg2 = get_cg(x_le_MAC,config2_class2)   

