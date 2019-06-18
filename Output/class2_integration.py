# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:43:30 2019

@author: Lisa
"""

'CLASS 2 WEIGHT ESTIMATION'
import inputs.concept_1 as c1
import modules.weight_class2.class2_weight as c2m
import modules.CG.class2_CG  as c2cg
import Output.connection_departments as detailedsizing

from inputs.constants import n_max,V_dive

'CLASS 2 WEIGHT EXECUTION AND CG LOCATION ESTIMATION'
#configuration 1 
config1_class2     = c2m.Class2_weight(1,c1.N_pax[0],c1.MTOW[0],max(c1.MTOW),c1.M_carried_canard_MZF[0],min(c1.M_MZF), n_max[0],V_dive[0],c1.M_fuel[0], max(c1.T_req), c1.l_f[0],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[0], detailedsizing.l_h, c1.S, 0., c1.b, 0., detailedsizing.S_v[0],detailedsizing.S_h[0],c1.Cr_t,0,c1.lambda_2_rad,detailedsizing.lambda_h_2_rad[0], detailedsizing.lambda_v_2_rad[0],0, c1.S_fus[0])     


config1_M_structural            =config1_class2.structural_mass()
config1_M_powerplant            =config1_class2.powerplant_mass()
config1_M_fixedeq               =config1_class2.fixed_equipment_mass()

config1_M_winggroup             =config1_class2.get_wing_group_mass()
config1_M_fuselagegroup         =config1_class2.get_fuselage_group_mass()


config1_class2_OEW              =config1_class2.OEW(config1_M_structural,config1_M_powerplant,config1_M_fixedeq)
#cg locations
config1_cg = c2cg.get_cg(detailedsizing.x_le_MAC,config1_class2)   
config1_cg_x=config1_cg.calc_x_cg()
config1_cg_y=config1_cg.calc_y_cg()
config1_cg_z=config1_cg.calc_z_cg()


#configuration 3
#masses
config3_class2     = c2m.Class2_weight(3,c1.N_pax[2],c1.MTOW[2],max(c1.MTOW),c1.M_carried_canard_MZF[2],min(c1.M_MZF), n_max[2],V_dive[2],c1.M_fuel[2], max(c1.T_req), c1.l_f[2],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[2], detailedsizing.l_h, c1.S, detailedsizing.S_c3, c1.b, detailedsizing.b_c3, detailedsizing.S_v, detailedsizing.S_h,c1.Cr_t, detailedsizing.Cr_t_c3, c1.lambda_2_rad,detailedsizing.lambda_h_2_rad, detailedsizing.lambda_v_2_rad,detailedsizing.lambda_c_le_rad3, c1.S_fus[2])

config3_M_structural            =config3_class2.structural_mass()
config3_M_powerplant            =config3_class2.powerplant_mass()
config3_M_fixedeq               =config3_class2.fixed_equipment_mass()

config3_M_winggroup             =config3_class2.get_wing_group_mass()
config3_M_fuselagegroup         =config3_class2.get_fuselage_group_mass()

config3_class2_OEW              =config3_class2.OEW(config3_M_structural,config3_M_powerplant,config3_M_fixedeq)

#cg locations
config3_cg = c2cg.get_cg(detailedsizing.x_le_MAC,config3_class2) 
config3_cg_x=config3_cg.calc_x_cg()
config3_cg_y=config3_cg.calc_y_cg()
config3_cg_z=config3_cg.calc_z_cg()


#configuration 2
config2_class2     = c2m.Class2_weight(2,c1.N_pax[1],c1.MTOW[1],max(c1.MTOW),c1.M_carried_canard_MZF[1],min(c1.M_MZF), n_max[1],V_dive[1],c1.M_fuel[1], max(c1.T_req), c1.l_f[1],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[1], detailedsizing.l_h, c1.S, detailedsizing.S_c2, c1.b, detailedsizing.b_c2, detailedsizing.S_v,detailedsizing.S_h,c1.Cr_t,detailedsizing.Cr_t_c2,c1.lambda_2_rad,detailedsizing.lambda_h_2_rad, detailedsizing.lambda_v_2_rad,detailedsizing.lambda_c_le_rad2, c1.S_fus[1])
config2_M_structural            =config2_class2.structural_mass()
config2_M_powerplant            =config2_class2.powerplant_mass()
config2_M_powerplant            =config3_M_powerplant
config2_M_fixedeq               =config3_M_fixedeq                                  #needed as the fixed equipment needs to be the same for both
   
config2_M_winggroup             =config2_class2.get_wing_group_mass()
config2_M_fuselagegroup         =config2_class2.get_fuselage_group_mass()

                                          
config2_class2_OEW              =config2_class2.OEW(config2_M_structural,config2_M_powerplant,config2_M_fixedeq)
#cg locations
config2_cg = c2cg.get_cg(detailedsizing.x_le_MAC,config2_class2)   
config2_cg_x=config2_cg.calc_x_cg()
config2_cg_y=config2_cg.calc_y_cg()
config2_cg_z=config2_cg.calc_z_cg()