# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


from inputs.concept_1 import x_le_MAC, l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout
from inputs.constants import *

from modules.Stability_lastresort.loaddiagram_detailed import Loading_diagram
from modules.main_class2 import config1_class2, config1_cg_x, config1_cg, config1_class2_OEW, config2_class2, config2_cg_x, config2_cg, config2_class2_OEW, config3_class2, config3_cg_x, config3_cg, config3_class2_OEW
from modules.CG.class2_CG import get_cg



x_le_MAC1_emp = [x_le_MAC[0] - 0.1 * l_f[0], x_le_MAC[1] - 0.1 * l_f[0], x_le_MAC[2] - 0.1 * l_f[0]]
x_le_MAC2_emp = [x_le_MAC[0] , x_le_MAC[1], x_le_MAC[2]]
x_le_MAC3_emp = [x_le_MAC[0] + 0.1 * l_f[0], x_le_MAC[1] + 0.1 * l_f[0], x_le_MAC[2] + 0.1 * l_f[0]]

"""Config 1"""
x_cg_config1_range_emp = [get_cg(x_le_MAC1_emp,config1_class2).calc_x_cg(),get_cg(x_le_MAC2_emp,config1_class2).calc_x_cg(),get_cg(x_le_MAC3_emp,config1_class2).calc_x_cg()]
x_cg_wing_config1_range_emp = [get_cg(x_le_MAC1_emp,config1_class2).x_cg_wing,get_cg(x_le_MAC2_emp,config1_class2).x_cg_wing,get_cg(x_le_MAC3_emp,config1_class2).x_cg_wing]

x_le_MAC_range_emp1 = [x_le_MAC1_emp[0], x_le_MAC2_emp[0], x_le_MAC3_emp[0]]
x_le_MAC_range_perc_emp1 = [x_le_MAC1_emp[0]/l_f[0], x_le_MAC2_emp[0]/l_f[0], x_le_MAC3_emp[0]/l_f[0]]

x_cg_config1_range_emp1 = [config1_cg_x - 0.1* l_f[0],config1_cg_x,config1_cg_x + 0.1* l_f[0]]
x_cg_wing_config1_range_emp1 = [config1_cg.x_cg_wing - 0.1* l_f[0], config1_cg.x_cg_wing, config1_cg.x_cg_wing + 0.1* l_f[0]]

"""Config 2"""
x_cg_config2_range_emp = [get_cg(x_le_MAC1_emp,config2_class2).calc_x_cg(),get_cg(x_le_MAC2_emp,config2_class2).calc_x_cg(),get_cg(x_le_MAC3_emp,config2_class2).calc_x_cg()]
x_cg_wing_config2_range_emp = [get_cg(x_le_MAC1_emp,config2_class2).x_cg_wing,get_cg(x_le_MAC2_emp,config2_class2).x_cg_wing,get_cg(x_le_MAC3_emp,config2_class2).x_cg_wing]

x_le_MAC_range_emp2 = [x_le_MAC1_emp[1], x_le_MAC2_emp[1], x_le_MAC3_emp[1]]
x_le_MAC_range_perc_emp2 = [x_le_MAC1_emp[1]/l_f[1], x_le_MAC2_emp[1]/l_f[1], x_le_MAC3_emp[1]/l_f[1]]

x_cg_config1_range_emp2 = [config1_cg_x - 0.1* l_f[0],config1_cg_x,config1_cg_x + 0.1* l_f[0]]
x_cg_wing_config1_range_emp2 = [config1_cg.x_cg_wing - 0.1* l_f[0], config1_cg.x_cg_wing, config1_cg.x_cg_wing + 0.1* l_f[0]]

"""Config 3"""
x_cg_config3_range_emp = [get_cg(x_le_MAC1_emp,config3_class2).calc_x_cg(),get_cg(x_le_MAC2_emp,config3_class2).calc_x_cg(),get_cg(x_le_MAC3_emp,config3_class2).calc_x_cg()]
x_cg_wing_config3_range_emp = [get_cg(x_le_MAC1_emp,config3_class2).x_cg_wing,get_cg(x_le_MAC2_emp,config3_class2).x_cg_wing,get_cg(x_le_MAC3_emp,config3_class2).x_cg_wing]

x_le_MAC_range_emp3 = [x_le_MAC1_emp[2], x_le_MAC2_emp[2], x_le_MAC3_emp[2]]
x_le_MAC_range_perc_emp3 = [x_le_MAC1_emp[2]/l_f[2], x_le_MAC2_emp[2]/l_f[2], x_le_MAC3_emp[2]/l_f[2]]

x_cg_config1_range_emp3 = [config3_cg_x - 0.1* l_f[0],config1_cg_x,config1_cg_x + 0.1* l_f[0]]
x_cg_wing_config1_range_emp3 = [config3_cg.x_cg_wing - 0.1* l_f[0], config1_cg.x_cg_wing, config1_cg.x_cg_wing + 0.1* l_f[0]]

"""Config 1"""
config1_load_emp       = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range_emp[0], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range_emp[0], -1, l_cutout)     
config1_load2_emp      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range_emp[1], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range_emp[1], 0, l_cutout)     
config1_load3_emp      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, x_cg_config1_range_emp[2], MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing_config1_range_emp[2], 1, l_cutout)     

"""Config 2"""
config2_load_emp       = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range_emp[0], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range_emp[0], -1, l_cutout)     
config2_load2_emp      = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range_emp[1], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range_emp[1], 0, l_cutout)     
config2_load3_emp      = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, x_cg_config2_range_emp[2], MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing_config2_range_emp[2], 1, l_cutout)     

"""Config 3"""
config3_load_emp       = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, config3_class2_OEW, x_cg_config3_range_emp[0], MAC, S, b, A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, x_cg_wing_config3_range_emp[0], -1, l_cutout)     
config3_load2_emp      = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, config3_class2_OEW, x_cg_config3_range_emp[1], MAC, S, b, A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, x_cg_wing_config3_range_emp[1], 0, l_cutout)     
config3_load3_emp      = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, config3_class2_OEW, x_cg_config3_range_emp[2], MAC, S, b, A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, x_cg_wing_config3_range_emp[2], 1, l_cutout)     


"""Config 1"""
cg1_pass_emp1 = [0, 0, 0]
cg2_pass_emp1 = [0, 0, 0]
weight_pass_emp1 = [0, 0, 0]
xcg_max_emp1 = [0, 0, 0, 0, 0, 0]
xcg_min_emp1 = [0, 0, 0, 0, 0, 0]
cg1_pass_emp1[0], cg2_pass_emp1[0], weight_pass_emp1[0], xcg_max_emp1[0], xcg_min_emp1[0] = config1_load_emp.loading_diagrams_pass()
cg1_pass_emp1[1], cg2_pass_emp1[1], weight_pass_emp1[1], xcg_max_emp1[1], xcg_min_emp1[1] = config1_load2_emp.loading_diagrams_pass()
cg1_pass_emp1[2], cg2_pass_emp1[2], weight_pass_emp1[2], xcg_max_emp1[2], xcg_min_emp1[2] = config1_load3_emp.loading_diagrams_pass()

cg1_fuel_emp1 = [0, 0, 0]
cg2_fuel_emp1 = [0, 0, 0]
weight_fuel_emp1 = [0, 0, 0]
cg1_fuel_emp1[0], cg2_fuel_emp1[0], weight_fuel_emp1[0], xcg_max_emp1[3], xcg_min_emp1[3] = config1_load_emp.loading_diagrams_fuel()
cg1_fuel_emp1[1], cg2_fuel_emp1[1], weight_fuel_emp1[1], xcg_max_emp1[4], xcg_min_emp1[4] = config1_load2_emp.loading_diagrams_fuel()
cg1_fuel_emp1[2], cg2_fuel_emp1[2], weight_fuel_emp1[2], xcg_max_emp1[5], xcg_min_emp1[5] = config1_load3_emp.loading_diagrams_fuel()

x_cg_max1_emp1 = [0, 0, 0]
x_cg_min1_emp1 = [0, 0, 0]

x_cg_max1_emp1[0] = max(xcg_max_emp1[0], xcg_max_emp1[3]) + 0.05*MAC
x_cg_max1_emp1[1] = max(xcg_max_emp1[1], xcg_max_emp1[4]) + 0.05*MAC
x_cg_max1_emp1[2] = max(xcg_max_emp1[2], xcg_max_emp1[5]) + 0.05*MAC
x_cg_min1_emp1[0] = min(xcg_min_emp1[0], xcg_min_emp1[3]) - 0.05*MAC
x_cg_min1_emp1[1] = min(xcg_min_emp1[1], xcg_min_emp1[4]) - 0.05*MAC
x_cg_min1_emp1[2] = min(xcg_min_emp1[2], xcg_min_emp1[5]) - 0.05*MAC

"""Config 2"""
cg1_pass_emp2 = [0, 0, 0]
cg2_pass_emp2 = [0, 0, 0]
weight_pass_emp2 = [0, 0, 0]
xcg_max_emp2 = [0, 0, 0, 0, 0, 0]
xcg_min_emp2 = [0, 0, 0, 0, 0, 0]
cg1_pass_emp2[0], cg2_pass_emp2[0], weight_pass_emp2[0], xcg_max_emp2[0], xcg_min_emp2[0] = config2_load_emp.loading_diagrams_pass()
cg1_pass_emp2[1], cg2_pass_emp2[1], weight_pass_emp2[1], xcg_max_emp2[1], xcg_min_emp2[1] = config2_load2_emp.loading_diagrams_pass()
cg1_pass_emp2[2], cg2_pass_emp2[2], weight_pass_emp2[2], xcg_max_emp2[2], xcg_min_emp2[2] = config2_load3_emp.loading_diagrams_pass()

cg1_fuel_emp2 = [0, 0, 0]
cg2_fuel_emp2 = [0, 0, 0]
weight_fuel_emp2 = [0, 0, 0]
cg1_fuel_emp2[0], cg2_fuel_emp2[0], weight_fuel_emp2[0], xcg_max_emp2[3], xcg_min_emp2[3] = config2_load_emp.loading_diagrams_fuel()
cg1_fuel_emp2[1], cg2_fuel_emp2[1], weight_fuel_emp2[1], xcg_max_emp2[4], xcg_min_emp2[4] = config2_load2_emp.loading_diagrams_fuel()
cg1_fuel_emp2[2], cg2_fuel_emp2[2], weight_fuel_emp2[2], xcg_max_emp2[5], xcg_min_emp2[5] = config2_load3_emp.loading_diagrams_fuel()

x_cg_max1_emp2 = [0, 0, 0]
x_cg_min1_emp2 = [0, 0, 0]

x_cg_max1_emp2[0] = max(xcg_max_emp2[0], xcg_max_emp2[3]) + 0.05*MAC
x_cg_max1_emp2[1] = max(xcg_max_emp2[1], xcg_max_emp2[4]) + 0.05*MAC
x_cg_max1_emp2[2] = max(xcg_max_emp2[2], xcg_max_emp2[5]) + 0.05*MAC
x_cg_min1_emp2[0] = min(xcg_min_emp2[0], xcg_min_emp2[3]) - 0.05*MAC
x_cg_min1_emp2[1] = min(xcg_min_emp2[1], xcg_min_emp2[4]) - 0.05*MAC
x_cg_min1_emp2[2] = min(xcg_min_emp2[2], xcg_min_emp2[5]) - 0.05*MAC

"""Config 3"""
cg1_pass_emp3 = [0, 0, 0]
cg2_pass_emp3 = [0, 0, 0]
weight_pass_emp3 = [0, 0, 0]
xcg_max_emp3 = [0, 0, 0, 0, 0, 0]
xcg_min_emp3 = [0, 0, 0, 0, 0, 0]
cg1_pass_emp3[0], cg2_pass_emp3[0], weight_pass_emp3[0], xcg_max_emp3[0], xcg_min_emp3[0] = config3_load_emp.loading_diagrams_pass()
cg1_pass_emp3[1], cg2_pass_emp3[1], weight_pass_emp3[1], xcg_max_emp3[1], xcg_min_emp3[1] = config3_load2_emp.loading_diagrams_pass()
cg1_pass_emp3[2], cg2_pass_emp3[2], weight_pass_emp3[2], xcg_max_emp3[2], xcg_min_emp3[2] = config3_load3_emp.loading_diagrams_pass()

cg1_fuel_emp3 = [0, 0, 0]
cg2_fuel_emp3 = [0, 0, 0]
weight_fuel_emp3 = [0, 0, 0]
cg1_fuel_emp3[0], cg2_fuel_emp3[0], weight_fuel_emp3[0], xcg_max_emp3[3], xcg_min_emp3[3] = config3_load_emp.loading_diagrams_fuel()
cg1_fuel_emp3[1], cg2_fuel_emp3[1], weight_fuel_emp3[1], xcg_max_emp3[4], xcg_min_emp3[4] = config3_load2_emp.loading_diagrams_fuel()
cg1_fuel_emp3[2], cg2_fuel_emp3[2], weight_fuel_emp3[2], xcg_max_emp3[5], xcg_min_emp3[5] = config3_load3_emp.loading_diagrams_fuel()

x_cg_max1_emp3 = [0, 0, 0]
x_cg_min1_emp3 = [0, 0, 0]

x_cg_max1_emp3[0] = max(xcg_max_emp3[0], xcg_max_emp3[3]) + 0.05*MAC
x_cg_max1_emp3[1] = max(xcg_max_emp3[1], xcg_max_emp3[4]) + 0.05*MAC
x_cg_max1_emp3[2] = max(xcg_max_emp3[2], xcg_max_emp3[5]) + 0.05*MAC
x_cg_min1_emp3[0] = min(xcg_min_emp3[0], xcg_min_emp3[3]) - 0.05*MAC
x_cg_min1_emp3[1] = min(xcg_min_emp3[1], xcg_min_emp3[4]) - 0.05*MAC
x_cg_min1_emp3[2] = min(xcg_min_emp3[2], xcg_min_emp3[5]) - 0.05*MAC
