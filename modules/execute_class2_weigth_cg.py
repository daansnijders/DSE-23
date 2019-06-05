# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:58:24 2019

@author: Lisa
"""
from modules.weight_class2.class2_weight import *
from modules.CG.class2_CG import *
from inputs.concept_1 import *



#class 2 execution  in terms of weight and cg 

#configuration 1 
config1_class2     = Class2_weight(1,N_pax[0],MTOW[0],M_carried_canard_MZF[0],min(M_MZF), n_max[0],V_dive[0],M_fuel[0], max(T_req), l_f[0],d_f_inner,d_f_outer,l_cabin[0], l_h[0], S, S_c[0], b, b_c[0], S_v[0],S_h[0],Cr_t,Cr_t_c[0],lambda_2_rad,lambda_h_2_rad[0], lambda_v_2_rad[0],lambda_c_2_rad, S_fus[0])     


config1_M_structural            =config1_class2.structural_mass()
config1_M_powerplant            =config1_class2.powerplant_mass()
config1_M_fixedeq               =config1_class2.fixed_equipment_mass()

config1_M_winggroup             =config1_class2.get_wing_group_mass()
config1_M_fuselagegroup         =config1_class2.get_fuselage_group_mass()

config1_class2_OEW              =config1_class2.OEW(config1_M_structural,config1_M_powerplant,config1_M_fixedeq)
#cg locations
config1_cg = get_cg(x_le_MAC,config1_class2)   
config1_cg_x=config1_cg.calc_x_cg()
config1_cg_y=config1_cg.calc_y_cg()
config1_cg_z=config1_cg.calc_z_cg()

#configuration 3
#masses
config3_class2     = Class2_weight(3,N_pax[2],MTOW[2],M_carried_canard_MZF[2],min(M_MZF), n_max[2],V_dive[2],M_fuel[2], max(T_req), l_f[2],d_f_inner,d_f_outer,l_cabin[2], l_h[2], S, S_c[2], b, b_c[2], S_v[2],S_h[2],Cr_t,Cr_t_c[2],lambda_2_rad,lambda_h_2_rad[2], lambda_v_2_rad[2],lambda_c_2_rad, S_fus[2])

config3_M_structural            =config3_class2.structural_mass()
config3_M_powerplant            =config3_class2.powerplant_mass()
config3_M_fixedeq               =config3_class2.fixed_equipment_mass()

config3_M_winggroup             =config3_class2.get_wing_group_mass()
config3_M_fuselagegroup         =config3_class2.get_fuselage_group_mass()

config3_class2_OEW              =config3_class2.OEW(config3_M_structural,config3_M_powerplant,config3_M_fixedeq)
#cg locations
config3_cg = get_cg(x_le_MAC,config3_class2) 
config3_cg_x=config3_cg.calc_x_cg()
config3_cg_y=config3_cg.calc_y_cg()
config3_cg_z=config3_cg.calc_z_cg()

#configuration 2
config2_class2     = Class2_weight(2,N_pax[1],MTOW[1],M_carried_canard_MZF[1],min(M_MZF), n_max[1],V_dive[1],M_fuel[1], max(T_req), l_f[1],d_f_inner,d_f_outer,l_cabin[1], l_h[1], S, S_c[1], b, b_c[1], S_v[1],S_h[1],Cr_t,Cr_t_c[1],lambda_2_rad,lambda_h_2_rad[1], lambda_v_2_rad[1],lambda_c_2_rad, S_fus[1])
config2_M_structural            =config2_class2.structural_mass()
config2_M_powerplant            =config2_class2.powerplant_mass()
config2_M_fixedeq               =config3_M_fixedeq                                  #needed as the fixed equipment needs to be the same for both
   
config2_M_winggroup             =config2_class2.get_wing_group_mass()
config2_M_fuselagegroup         =config2_class2.get_fuselage_group_mass()
                                          
config2_class2_OEW              =config2_class2.OEW(config2_M_structural,config2_M_powerplant,config2_M_fixedeq)
#cg locations
config2_cg = get_cg(x_le_MAC,config2_class2)   
config2_cg_x=config2_cg.calc_x_cg()
config2_cg_y=config2_cg.calc_y_cg()
config2_cg_z=config2_cg.calc_z_cg()


