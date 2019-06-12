# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


from inputs.concept_1 import *
from inputs.constants import *

from modules.Stability.loaddiagram_detailed import *
from modules.Stability.Check_ground import *
from modules.main_class2 import *

x_le_MAC1 = [x_le_MAC[0] - 0.1 * l_f[0], x_le_MAC[1] - 0.1 * l_f[1], x_le_MAC[2] - 0.1 * l_f[2]]
x_le_MAC2 = [x_le_MAC[0] , x_le_MAC[1], x_le_MAC[2]]
x_le_MAC3 = [x_le_MAC[0] + 0.1 * l_f[0], x_le_MAC[1] + 0.1 * l_f[1], x_le_MAC[2] + 0.1 * l_f[2]]

x_cg_config1_range = [get_cg(x_le_MAC1,config1_class2).calc_x_cg(),get_cg(x_le_MAC2,config1_class2).calc_x_cg(),get_cg(x_le_MAC3,config1_class2).calc_x_cg()]
x_cg_wing_config1_range = [get_cg(x_le_MAC1,config1_class2).x_cg_wing,get_cg(x_le_MAC2,config1_class2).x_cg_wing,get_cg(x_le_MAC3,config1_class2).x_cg_wing]

x_le_MAC_range = [x_le_MAC1[0], x_le_MAC2[0], x_le_MAC3[0]]
x_le_MAC_range_perc = [x_le_MAC1[0]/l_f[0], x_le_MAC2[0]/l_f[0], x_le_MAC3[0]/l_f[0]]

#config1_cg_x=config1_cg.calc_x_cg()

x_cg_config1_range = [config1_cg_x - 0.1* l_f[0],config1_cg_x,config1_cg_x + 0.1* l_f[0]]
x_cg_wing_config1_range = [config1_cg.x_cg_wing - 0.1* l_f[0], config1_cg.x_cg_wing, config1_cg.x_cg_wing + 0.1* l_f[0]]


config1_load      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range[0], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range[0], -1)     
config1_load2      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range[1], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range[1], 0)     
config1_load3      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range[2], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range[2], 1)     


cg1_pass = [0, 0, 0]
cg2_pass = [0, 0, 0]
weight_pass = [0, 0, 0]
xcg_max = [0, 0, 0, 0, 0, 0]
xcg_min = [0, 0, 0, 0, 0, 0]
cg1_pass[0], cg2_pass[0], weight_pass[0], xcg_max[0], xcg_min[0] = config1_load.loading_diagrams_pass()
cg1_pass[1], cg2_pass[1], weight_pass[1], xcg_max[1], xcg_min[1] = config1_load2.loading_diagrams_pass()
cg1_pass[2], cg2_pass[2], weight_pass[2], xcg_max[2], xcg_min[2] = config1_load3.loading_diagrams_pass()

cg1_fuel = [0, 0, 0]
cg2_fuel = [0, 0, 0]
weight_fuel = [0, 0, 0]
cg1_fuel[0], cg2_fuel[0], weight_fuel[0], xcg_max[3], xcg_min[3] = config1_load.loading_diagrams_fuel()
cg1_fuel[1], cg2_fuel[1], weight_fuel[1], xcg_max[4], xcg_min[4] = config1_load2.loading_diagrams_fuel()
cg1_fuel[2], cg2_fuel[2], weight_fuel[2], xcg_max[5], xcg_min[5] = config1_load3.loading_diagrams_fuel()


#config1_ground       = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config1_ground2      = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config1_ground3      = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     


#frac_min = [0,0,0]
#frac_max = [0,0,0]
#frac_min[0], frac_max[0], frac1 = config1_ground.check_equilibrium()
#frac_min[1], frac_max[1], frac2 = config1_ground2.check_equilibrium()
#frac_min[2], frac_max[2], frac2 = config1_ground3.check_equilibrium()

x_cg_max1 = [0, 0, 0]
x_cg_min1 = [0, 0, 0]

x_cg_max1[0] = max(xcg_max[0], xcg_max[3]) + 0.05*MAC
x_cg_max1[1] = max(xcg_max[1], xcg_max[4]) + 0.05*MAC
x_cg_max1[2] = max(xcg_max[2], xcg_max[5]) + 0.05*MAC
x_cg_min1[0] = min(xcg_min[0], xcg_min[3]) - 0.05*MAC
x_cg_min1[1] = min(xcg_min[1], xcg_min[4]) - 0.05*MAC
x_cg_min1[2] = min(xcg_min[2], xcg_min[5]) - 0.05*MAC

#print ("The most aft CG position from the nose for configuration 1 during flight is: ", max(xcg_max[0], xcg_max[3]))
#print("Including a 0,05 m stability margin we get", max(xcg_max)+0.05)
#
#print ("The most forward CG position from the nose for configuration 1 during flight is: ", min(xcg_min[0], xcg_min[3]))



