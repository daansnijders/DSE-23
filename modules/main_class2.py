# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:58:24 2019

@author: Lisa
"""

import inputs.concept_1 as c1
import modules.weight_class2.class2_weight as c2m
import modules.CG.class2_CG  as c2cg

from inputs.constants import n_max,V_dive

'CLASS 2 WEIGHT EXECUTION AND CG LOCATION ESTIMATION'
#configuration 1 
config1_class2     = c2m.Class2_weight(1,c1.N_pax[0],c1.MTOW[0],max(c1.MTOW),c1.M_carried_canard_MZF[0],min(c1.M_MZF), n_max[0],V_dive[0],c1.M_fuel[0], max(c1.T_req), c1.l_f[0],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[0], c1.l_h[0], c1.S, c1.S_c[0], c1.b, c1.b_c[0], c1.S_v[0],c1.S_h[0],c1.Cr_t,c1.Cr_t_c[0],c1.lambda_2_rad,c1.lambda_h_2_rad[0], c1.lambda_v_2_rad[0],c1.lambda_c_2_rad, c1.S_fus[0])     


config1_M_structural            =config1_class2.structural_mass()
config1_M_powerplant            =config1_class2.powerplant_mass()
config1_M_fixedeq               =config1_class2.fixed_equipment_mass()

config1_M_winggroup             =config1_class2.get_wing_group_mass()
config1_M_fuselagegroup         =config1_class2.get_fuselage_group_mass()
#x_le_MAC=get_x_le_MAC(l_f,MAC,M_wing_group, M_fuselage_group)

config1_class2_OEW              =config1_class2.OEW(config1_M_structural,config1_M_powerplant,config1_M_fixedeq)
#cg locations
config1_cg = c2cg.get_cg(c1.x_le_MAC,config1_class2,c1.b_h[0],c1.Cr_h[0],c1.Ct_h[0],c1.lambda_h_le_rad,c1.x_le_h[0],c1.b_v[0],c1.Cr_v[0],c1.Ct_v[0],c1.lambda_v_le_rad,c1.x_le_v[0],c1.Cr_c[0],c1.t_c_c[0],c1.z_mlg,c1.z_nlg,c1.x_mlg,c1.x_nlg)   
config1_cg_x=config1_cg.calc_x_cg()
config1_cg_y=config1_cg.calc_y_cg()
config1_cg_z=config1_cg.calc_z_cg()


#configuration 3
#masses
config3_class2     = c2m.Class2_weight(3,c1.N_pax[2],c1.MTOW[2],max(c1.MTOW),c1.M_carried_canard_MZF[2],min(c1.M_MZF), n_max[2],V_dive[2],c1.M_fuel[2], max(c1.T_req), c1.l_f[2],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[2], c1.l_h[2], c1.S, c1.S_c[2], c1.b, c1.b_c[2], c1.S_v[2],c1.S_h[2],c1.Cr_t,c1.Cr_t_c[2],c1.lambda_2_rad,c1.lambda_h_2_rad[2], c1.lambda_v_2_rad[2],c1.lambda_c_2_rad, c1.S_fus[2])

config3_M_structural            =config3_class2.structural_mass()
config3_M_powerplant            =config3_class2.powerplant_mass()
config3_M_fixedeq               =config3_class2.fixed_equipment_mass()

config3_M_winggroup             =config3_class2.get_wing_group_mass()
config3_M_fuselagegroup         =config3_class2.get_fuselage_group_mass()

config3_class2_OEW              =config3_class2.OEW(config3_M_structural,config3_M_powerplant,config3_M_fixedeq)

#cg locations
config3_cg = c2cg.get_cg(c1.x_le_MAC,config3_class2,c1.b_h[2],c1.Cr_h[2],c1.Ct_h[2],c1.lambda_h_le_rad,c1.x_le_h[2],c1.b_v[2],c1.Cr_v[2],c1.Ct_v[2],c1.lambda_v_le_rad,c1.x_le_v[2],c1.Cr_c[2],c1.t_c_c[2],c1.z_mlg,c1.z_nlg,c1.x_mlg,c1.x_nlg) 
config3_cg_x=config3_cg.calc_x_cg()
config3_cg_y=config3_cg.calc_y_cg()
config3_cg_z=config3_cg.calc_z_cg()


#configuration 2
config2_class2     = c2m.Class2_weight(2,c1.N_pax[1],c1.MTOW[1],max(c1.MTOW),c1.M_carried_canard_MZF[1],min(c1.M_MZF), n_max[1],V_dive[1],c1.M_fuel[1], max(c1.T_req), c1.l_f[1],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[1], c1.l_h[1], c1.S, c1.S_c[1], c1.b, c1.b_c[1], c1.S_v[1],c1.S_h[1],c1.Cr_t,c1.Cr_t_c[1],c1.lambda_2_rad,c1.lambda_h_2_rad[1], c1.lambda_v_2_rad[1],c1.lambda_c_2_rad, c1.S_fus[1])
config2_M_structural            =config2_class2.structural_mass()
config2_M_powerplant            =config2_class2.powerplant_mass()
config2_M_powerplant            =config3_M_powerplant
config2_M_fixedeq               =config3_M_fixedeq                                  #needed as the fixed equipment needs to be the same for both
   
config2_M_winggroup             =config2_class2.get_wing_group_mass()
config2_M_fuselagegroup         =config2_class2.get_fuselage_group_mass()

                                          
config2_class2_OEW              =config2_class2.OEW(config2_M_structural,config2_M_powerplant,config2_M_fixedeq)
#cg locations
config2_cg = c2cg.get_cg(c1.x_le_MAC,config2_class2,c1.b_h[1],c1.Cr_h[1],c1.Ct_h[1],c1.lambda_h_le_rad,c1.x_le_h[1],c1.b_v[1],c1.Cr_v[1],c1.Ct_v[1],c1.lambda_v_le_rad,c1.x_le_v[1],c1.Cr_c[1],c1.t_c_c[1],c1.z_mlg,c1.z_nlg,c1.x_mlg,c1.x_nlg) 
 
config2_cg_x=config2_cg.calc_x_cg()
config2_cg_y=config2_cg.calc_y_cg()
config2_cg_z=config2_cg.calc_z_cg()





