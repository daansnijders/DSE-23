# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


#from inputs.concept_1 import x_le_MAC, l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout
from inputs.concept_1 import l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout, x_le_h, x_le_v #new iteration
from inputs.constants import *
import Output.read_load_variables as varib

from modules.CG.class2_CG import get_cg
from modules.Stability.loaddiagram_detailed import Loading_diagram
from Output.class2_integration import config3_class2, config3_cg_x, config3_cg, config3_class2_OEW

x_le_MAC = varib.x_le_MAC
x_le_MAC1_can2 = [x_le_MAC[0] - 0.1 * l_f[0], x_le_MAC[1] - 0.1 * l_f[1], x_le_MAC[2] - 0.1 * l_f[2]]
x_le_MAC2_can2 = [x_le_MAC[0] , x_le_MAC[1], x_le_MAC[2]]
x_le_MAC3_can2 = [x_le_MAC[0] + 0.1 * l_f[0], x_le_MAC[1] + 0.1 * l_f[1], x_le_MAC[2] + 0.1 * l_f[2]]


config3_cgrange1=get_cg(x_le_MAC1_can2,config3_class2,varib.b_h,varib.Cr_h,varib.Ct_h,varib.lambda_h_le_rad,x_le_h[2],varib.b_v,varib.Cr_v,varib.Ct_v,varib.lambda_v_le_rad,x_le_v[2],varib.Cr_c3,varib.t_c_c3,varib.z_mlg,varib.z_nlg,varib.x_mlg,varib.x_nlg)
config3_cgrange2=get_cg(x_le_MAC2_can2,config3_class2,varib.b_h,varib.Cr_h,varib.Ct_h,varib.lambda_h_le_rad,x_le_h[2],varib.b_v,varib.Cr_v,varib.Ct_v,varib.lambda_v_le_rad,x_le_v[2],varib.Cr_c3,varib.t_c_c3,varib.z_mlg,varib.z_nlg,varib.x_mlg,varib.x_nlg)
config3_cgrange3=get_cg(x_le_MAC3_can2,config3_class2,varib.b_h,varib.Cr_h,varib.Ct_h,varib.lambda_h_le_rad,x_le_h[2],varib.b_v,varib.Cr_v,varib.Ct_v,varib.lambda_v_le_rad,x_le_v[2],varib.Cr_c3,varib.t_c_c3,varib.z_mlg,varib.z_nlg,varib.x_mlg,varib.x_nlg)

x_cg_config3_range_can2 = [config3_cgrange1.calc_x_cg(),config3_cgrange2.calc_x_cg(),config3_cgrange3.calc_x_cg()]
x_cg_wing_config3_range_can2 = [config3_cgrange1.x_cg_wing,config3_cgrange2.x_cg_wing,config3_cgrange3.x_cg_wing]


#x_cg_config3_range_can2 = [get_cg(x_le_MAC1_can2,config3_class2).calc_x_cg(),get_cg(x_le_MAC2_can2,config3_class2).calc_x_cg(),get_cg(x_le_MAC3_can2,config3_class2).calc_x_cg()]
#x_cg_wing_config3_range_can2 = [get_cg(x_le_MAC1_can2,config3_class2).x_cg_wing,get_cg(x_le_MAC2_can2,config3_class2).x_cg_wing,get_cg(x_le_MAC3_can2,config3_class2).x_cg_wing]

x_le_MAC_range_can2 = [x_le_MAC1_can2[1], x_le_MAC2_can2[1], x_le_MAC3_can2[1]]
#x_le_MAC_range_perc = [x_le_MAC1[1]/l_f[1], x_le_MAC2[1]/l_f[1], x_le_MAC3[1]/l_f[1]]
x_le_MAC_range_perccanard3_can2 = [x_le_MAC2_can2[1]/l_f[1]]

#config2_cg_x=config2_cg.calc_x_cg()

x_cg_config2_range_can2 = [config3_cg_x - 0.1* l_f[1],config3_cg_x,config3_cg_x + 0.1* l_f[1]]
x_cg_wing_config2_range_can2 = [config3_cg.x_cg_wing - 0.1* l_f[1], config3_cg.x_cg_wing, config3_cg.x_cg_wing + 0.1* l_f[1]]


#config2_load       = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range[0], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range[0], -1, l_cutout)     
config2_load2_can2      = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, config3_class2_OEW, x_cg_config3_range_can2[1], MAC, S, b, A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, x_cg_wing_config3_range_can2[1], 0, l_cutout)     
#config2_load3      = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range[2], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range[2], 1, l_cutout)     


cg1_pass_can2 = [0, 0, 0]
cg2_pass_can2 = [0, 0, 0]
weight_pass_can2 = [0, 0, 0]
xcg_max_can2 = [0, 0, 0, 0, 0, 0]
xcg_min_can2 = [0, 0, 0, 0, 0, 0]
#cg1_pass[0], cg2_pass[0], weight_pass[0], xcg_max[0], xcg_min[0] = config2_load.loading_diagrams_pass()
cg1_pass_can2[1], cg2_pass_can2[1], weight_pass_can2[1], xcg_max_can2[1], xcg_min_can2[1] = config2_load2_can2.loading_diagrams_pass()
#cg1_pass[2], cg2_pass[2], weight_pass[2], xcg_max[2], xcg_min[2] = config2_load3.loading_diagrams_pass()

cg1_fuel_can2 = [0, 0, 0]
cg2_fuel_can2 = [0, 0, 0]
weight_fuel_can2 = [0, 0, 0]
#cg1_fuel[0], cg2_fuel[0], weight_fuel[0], xcg_max[3], xcg_min[3] = config2_load.loading_diagrams_fuel()
cg1_fuel_can2[1], cg2_fuel_can2[1], weight_fuel_can2[1], xcg_max_can2[4], xcg_min_can2[4] = config2_load2_can2.loading_diagrams_fuel()
#cg1_fuel[2], cg2_fuel[2], weight_fuel[2], xcg_max[5], xcg_min[5] = config2_load3.loading_diagrams_fuel()


#config2_ground       = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config2_ground2      = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config2_ground3      = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     


#frac_min = [0,0,0]
#frac_max = [0,0,0]
#frac_min[0], frac_max[0], frac1 = config2_ground.check_equilibrium()
#frac_min[1], frac_max[1], frac2 = config2_ground2.check_equilibrium()
#frac_min[2], frac_max[2], frac2 = config2_ground3.check_equilibrium()

x_cg_max3canard_can2 = [0]
x_cg_min3canard_can2 = [0]

#x_cg_max1canard[0] = max(xcg_max[0], xcg_max[3]) + 0.05*MAC
x_cg_max3canard_can2[0] = max(xcg_max_can2[1], xcg_max_can2[4]) + 0.05*MAC
#x_cg_max1canard[2] = max(xcg_max[2], xcg_max[5]) + 0.05*MAC
#x_cg_min1canard[0] = min(xcg_min[0], xcg_min[3]) - 0.05*MAC
x_cg_min3canard_can2[0] = min(xcg_min_can2[1], xcg_min_can2[4]) - 0.05*MAC
#x_cg_min1canard[2] = min(xcg_min[2], xcg_min[5]) - 0.05*MAC

#print ("The most aft CG position from the nose for configuration 1 during flight is: ", max(xcg_max[0], xcg_max[3]))
#print("Including a 0,05 m stability margin we get", max(xcg_max)+0.05)
#
#print ("The most forward CG position from the nose for configuration 1 during flight is: ", min(xcg_min[0], xcg_min[3]))




