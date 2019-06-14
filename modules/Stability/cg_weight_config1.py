# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


from inputs.concept_1 import x_le_MAC, l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout
from inputs.constants import *

from modules.Stability.loaddiagram_detailed import Loading_diagram
from modules.main_class2 import get_cg, config1_class2, config1_cg_x, config1_cg, config1_class2_OEW

x_le_MAC1_emp = [x_le_MAC[0] - 0.1 * l_f[0], x_le_MAC[1] - 0.1 * l_f[1], x_le_MAC[2] - 0.1 * l_f[2]]
x_le_MAC2_emp = [x_le_MAC[0] , x_le_MAC[1], x_le_MAC[2]]
x_le_MAC3_emp = [x_le_MAC[0] + 0.1 * l_f[0], x_le_MAC[1] + 0.1 * l_f[1], x_le_MAC[2] + 0.1 * l_f[2]]

x_cg_config1_range_emp = [get_cg(x_le_MAC1_emp,config1_class2).calc_x_cg(),get_cg(x_le_MAC2_emp,config1_class2).calc_x_cg(),get_cg(x_le_MAC3_emp,config1_class2).calc_x_cg()]
x_cg_wing_config1_range_emp = [get_cg(x_le_MAC1_emp,config1_class2).x_cg_wing,get_cg(x_le_MAC2_emp,config1_class2).x_cg_wing,get_cg(x_le_MAC3_emp,config1_class2).x_cg_wing]

x_le_MAC_range_emp = [x_le_MAC1_emp[0], x_le_MAC2_emp[0], x_le_MAC3_emp[0]]
x_le_MAC_range_perc_emp = [x_le_MAC1_emp[0]/l_f[0], x_le_MAC2_emp[0]/l_f[0], x_le_MAC3_emp[0]/l_f[0]]

#config1_cg_x=config1_cg.calc_x_cg()

x_cg_config1_range_emp = [config1_cg_x - 0.1* l_f[0],config1_cg_x,config1_cg_x + 0.1* l_f[0]]
x_cg_wing_config1_range_emp = [config1_cg.x_cg_wing - 0.1* l_f[0], config1_cg.x_cg_wing, config1_cg.x_cg_wing + 0.1* l_f[0]]


config1_load_emp      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range_emp[0], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range_emp[0], -1, l_cutout)     
config1_load2_emp      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range_emp[1], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range_emp[1], 0, l_cutout)     
config1_load3_emp      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range_emp[2], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range_emp[2], 1, l_cutout)     


cg1_pass_emp = [0, 0, 0]
cg2_pass_emp = [0, 0, 0]
weight_pass_emp = [0, 0, 0]
xcg_max_emp = [0, 0, 0, 0, 0, 0]
xcg_min_emp = [0, 0, 0, 0, 0, 0]
cg1_pass_emp[0], cg2_pass_emp[0], weight_pass_emp[0], xcg_max_emp[0], xcg_min_emp[0] = config1_load_emp.loading_diagrams_pass()
cg1_pass_emp[1], cg2_pass_emp[1], weight_pass_emp[1], xcg_max_emp[1], xcg_min_emp[1] = config1_load2_emp.loading_diagrams_pass()
cg1_pass_emp[2], cg2_pass_emp[2], weight_pass_emp[2], xcg_max_emp[2], xcg_min_emp[2] = config1_load3_emp.loading_diagrams_pass()

cg1_fuel_emp = [0, 0, 0]
cg2_fuel_emp = [0, 0, 0]
weight_fuel_emp = [0, 0, 0]
cg1_fuel_emp[0], cg2_fuel_emp[0], weight_fuel_emp[0], xcg_max_emp[3], xcg_min_emp[3] = config1_load_emp.loading_diagrams_fuel()
cg1_fuel_emp[1], cg2_fuel_emp[1], weight_fuel_emp[1], xcg_max_emp[4], xcg_min_emp[4] = config1_load2_emp.loading_diagrams_fuel()
cg1_fuel_emp[2], cg2_fuel_emp[2], weight_fuel_emp[2], xcg_max_emp[5], xcg_min_emp[5] = config1_load3_emp.loading_diagrams_fuel()

x_cg_max1_emp = [0, 0, 0]
x_cg_min1_emp = [0, 0, 0]

x_cg_max1_emp[0] = max(xcg_max_emp[0], xcg_max_emp[3]) + 0.05*MAC
x_cg_max1_emp[1] = max(xcg_max_emp[1], xcg_max_emp[4]) + 0.05*MAC
x_cg_max1_emp[2] = max(xcg_max_emp[2], xcg_max_emp[5]) + 0.05*MAC
x_cg_min1_emp[0] = min(xcg_min_emp[0], xcg_min_emp[3]) - 0.05*MAC
x_cg_min1_emp[1] = min(xcg_min_emp[1], xcg_min_emp[4]) - 0.05*MAC
x_cg_min1_emp[2] = min(xcg_min_emp[2], xcg_min_emp[5]) - 0.05*MAC


