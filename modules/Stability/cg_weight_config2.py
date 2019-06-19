# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


from inputs.concept_1 import x_le_MAC, l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout
#from inputs.concept_1 import l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout #new iteration
from inputs.constants import *

from modules.Stability.loaddiagram_detailed import Loading_diagram
from modules.main_class2 import config2_class2, config2_cg_x, config2_cg, config2_class2_OEW
from modules.CG.class2_CG import get_cg

x_le_MAC1_can1 = [x_le_MAC[0] - 0.1 * l_f[0], x_le_MAC[1] - 0.1 * l_f[1], x_le_MAC[2] - 0.1 * l_f[2]]
x_le_MAC2_can1 = [x_le_MAC[0] , x_le_MAC[1], x_le_MAC[2]]
x_le_MAC3_can1 = [x_le_MAC[0] + 0.1 * l_f[0], x_le_MAC[1] + 0.1 * l_f[1], x_le_MAC[2] + 0.1 * l_f[2]]

x_cg_config2_range_can1 = [get_cg(x_le_MAC1_can1,config2_class2).calc_x_cg(),get_cg(x_le_MAC2_can1,config2_class2).calc_x_cg(),get_cg(x_le_MAC3_can1,config2_class2).calc_x_cg()]
x_cg_wing_config2_range_can1 = [get_cg(x_le_MAC1_can1,config2_class2).x_cg_wing,get_cg(x_le_MAC2_can1,config2_class2).x_cg_wing,get_cg(x_le_MAC3_can1,config2_class2).x_cg_wing]

x_le_MAC_range_can1 = [x_le_MAC1_can1[1], x_le_MAC2_can1[1], x_le_MAC3_can1[1]]
#x_le_MAC_range_perc = [x_le_MAC1[1]/l_f[1], x_le_MAC2[1]/l_f[1], x_le_MAC3[1]/l_f[1]]
x_le_MAC_range_perccanard2_can1 = [x_le_MAC2_can1[1]/l_f[1]]

#config2_cg_x=config2_cg.calc_x_cg()

x_cg_config2_range_can1 = [config2_cg_x - 0.1* l_f[1],config2_cg_x,config2_cg_x + 0.1* l_f[1]]
x_cg_wing_config2_range_can1 = [config2_cg.x_cg_wing - 0.1* l_f[1], config2_cg.x_cg_wing, config2_cg.x_cg_wing + 0.1* l_f[1]]


#config2_load       = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range[0], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range[0], -1, l_cutout)     
config2_load2_can1      = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range_can1[1], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range_can1[1], 0, l_cutout)     
#config2_load3      = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range[2], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range[2], 1, l_cutout)     


cg1_pass_can1 = [0, 0, 0]
cg2_pass_can1 = [0, 0, 0]
weight_pass_can1 = [0, 0, 0]
xcg_max_can1 = [0, 0, 0, 0, 0, 0]
xcg_min_can1 = [0, 0, 0, 0, 0, 0]
#cg1_pass[0], cg2_pass[0], weight_pass[0], xcg_max[0], xcg_min[0] = config2_load.loading_diagrams_pass()
cg1_pass_can1[1], cg2_pass_can1[1], weight_pass_can1[1], xcg_max_can1[1], xcg_min_can1[1] = config2_load2_can1.loading_diagrams_pass()
#cg1_pass[2], cg2_pass[2], weight_pass[2], xcg_max[2], xcg_min[2] = config2_load3.loading_diagrams_pass()

cg1_fuel_can1 = [0, 0, 0]
cg2_fuel_can1 = [0, 0, 0]
weight_fuel_can1 = [0, 0, 0]
#cg1_fuel[0], cg2_fuel[0], weight_fuel[0], xcg_max[3], xcg_min[3] = config2_load.loading_diagrams_fuel()
cg1_fuel_can1[1], cg2_fuel_can1[1], weight_fuel_can1[1], xcg_max_can1[4], xcg_min_can1[4] = config2_load2_can1.loading_diagrams_fuel()
#cg1_fuel[2], cg2_fuel[2], weight_fuel[2], xcg_max[5], xcg_min[5] = config2_load3.loading_diagrams_fuel()


#config2_ground       = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config2_ground2      = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config2_ground3      = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     


#frac_min = [0,0,0]
#frac_max = [0,0,0]
#frac_min[0], frac_max[0], frac1 = config2_ground.check_equilibrium()
#frac_min[1], frac_max[1], frac2 = config2_ground2.check_equilibrium()
#frac_min[2], frac_max[2], frac2 = config2_ground3.check_equilibrium()

x_cg_max2canard_can1 = [0]
x_cg_min2canard_can1 = [0]

#x_cg_max1canard[0] = max(xcg_max[0], xcg_max[3]) + 0.05*MAC
x_cg_max2canard_can1[0] = max(xcg_max_can1[1], xcg_max_can1[4]) + 0.05*MAC
#x_cg_max1canard[2] = max(xcg_max[2], xcg_max[5]) + 0.05*MAC
#x_cg_min1canard[0] = min(xcg_min[0], xcg_min[3]) - 0.05*MAC
x_cg_min2canard_can1[0] = min(xcg_min_can1[1], xcg_min_can1[4]) - 0.05*MAC
#x_cg_min1canard[2] = min(xcg_min[2], xcg_min[5]) - 0.05*MAC

#print ("The most aft CG position from the nose for configuration 1 during flight is: ", max(xcg_max[0], xcg_max[3]))
#print("Including a 0,05 m stability margin we get", max(xcg_max)+0.05)
#
#print ("The most forward CG position from the nose for configuration 1 during flight is: ", min(xcg_min[0], xcg_min[3]))




