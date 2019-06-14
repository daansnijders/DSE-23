# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""

from inputs.concept_1 import *
from inputs.constants import *

from modules.Stability.loaddiagram_detailed import *
from modules.main_class2 import *

x_cg_config1_range = [config1_cg_x - 0.1* l_f[0],config1_cg_x,config1_cg_x + 0.1* l_f[0]]
x_cg_wing_config1_range = [config1_cg.x_cg_wing - 0.1* l_f[0], config1_cg.x_cg_wing, config1_cg.x_cg_wing + 0.1* l_f[0]]


config1_load      = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, config1_class2_OEW, config1_cg_x, MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, config1_cg.x_cg_wing, 1, l_cutout)     
config2_load      = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, config2_class2_OEW, config2_cg_x, MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, config2_cg.x_cg_wing, 2, l_cutout)     
config3_load      = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, config3_class2_OEW, config3_cg_x, MAC, S, b, A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, config3_cg.x_cg_wing, 3, l_cutout)     

cg1_pass = [0, 0, 0]
cg2_pass = [0, 0, 0]
weight_pass = [0, 0, 0]
xcg_max = [0, 0, 0, 0, 0, 0]
xcg_min = [0, 0, 0, 0, 0, 0]
cg1_pass[0], cg2_pass[0], weight_pass[0], xcg_max[0], xcg_min[0] = config1_load.loading_diagrams_pass()
cg1_pass[1], cg2_pass[1], weight_pass[1], xcg_max[1], xcg_min[1] = config2_load.loading_diagrams_pass()
cg1_pass[2], cg2_pass[2], weight_pass[2], xcg_max[2], xcg_min[2] = config3_load.loading_diagrams_pass()

cg1_fuel = [0, 0, 0]
cg2_fuel = [0, 0, 0]
weight_fuel = [0, 0, 0]
cg1_fuel[0], cg2_fuel[0], weight_fuel[0], xcg_max[3], xcg_min[3] = config1_load.loading_diagrams_fuel()
cg1_fuel[1], cg2_fuel[1], weight_fuel[1], xcg_max[4], xcg_min[4] = config2_load.loading_diagrams_fuel()
cg1_fuel[2], cg2_fuel[2], weight_fuel[2], xcg_max[5], xcg_min[5] = config3_load.loading_diagrams_fuel()

x_cg_max_flight1 = max(xcg_max[0], xcg_max[3]) + 0.05*MAC
x_cg_max_flight2 = max(xcg_max[1], xcg_max[4]) + 0.05*MAC
x_cg_max_flight3 = max(xcg_max[2], xcg_max[5]) + 0.05*MAC
x_cg_min_flight1 = min(xcg_min[0], xcg_min[3]) - 0.05*MAC
x_cg_min_flight2 = min(xcg_min[1], xcg_min[4]) - 0.05*MAC
x_cg_min_flight3 = min(xcg_min[2], xcg_min[5]) - 0.05*MAC





