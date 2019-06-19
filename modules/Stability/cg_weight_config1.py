# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


#from inputs.concept_1 import x_le_MAC, l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout

from inputs.concept_1 import l_f, x_cargo, l_cabin, N_pax, MAC, S, b, A, M_payload, M_cargo_available, M_fuel, l_cutout, x_le_h, x_le_v #new iteration
from inputs.constants import *
import Output.read_load_variables as varib

x_le_MAC = varib.x_le_MAC

from modules.Stability.loaddiagram_detailed import Loading_diagram
from Output.class2_integration import config1_class2, config1_cg_x, config1_cg, config1_class2_OEW
from modules.CG.class2_CG import get_cg

x_le_MAC1_emp = [x_le_MAC[0] - 0.1 * l_f[0], x_le_MAC[1] - 0.1 * l_f[1], x_le_MAC[2] - 0.1 * l_f[2]]
x_le_MAC2_emp = [x_le_MAC[0] , x_le_MAC[1], x_le_MAC[2]]
x_le_MAC3_emp = [x_le_MAC[0] + 0.1 * l_f[0], x_le_MAC[1] + 0.1 * l_f[1], x_le_MAC[2] + 0.1 * l_f[2]]


config1_cgrange1=get_cg(x_le_MAC1_emp,config1_class2,varib.b_h,varib.Cr_h,varib.Ct_h,varib.lambda_h_le_rad,x_le_h[0],varib.b_v,varib.Cr_v,varib.Ct_v,varib.lambda_v_le_rad,x_le_v[0],0,0,varib.z_mlg,varib.z_nlg,varib.x_mlg,varib.x_nlg)
config1_cgrange2=get_cg(x_le_MAC2_emp,config1_class2,varib.b_h,varib.Cr_h,varib.Ct_h,varib.lambda_h_le_rad,x_le_h[0],varib.b_v,varib.Cr_v,varib.Ct_v,varib.lambda_v_le_rad,x_le_v[0],0,0,varib.z_mlg,varib.z_nlg,varib.x_mlg,varib.x_nlg)
config1_cgrange3=get_cg(x_le_MAC3_emp,config1_class2,varib.b_h,varib.Cr_h,varib.Ct_h,varib.lambda_h_le_rad,x_le_h[0],varib.b_v,varib.Cr_v,varib.Ct_v,varib.lambda_v_le_rad,x_le_v[0],0,0,varib.z_mlg,varib.z_nlg,varib.x_mlg,varib.x_nlg)


x_cg_config1_range_emp = [config1_cgrange1.calc_x_cg(),config1_cgrange2.calc_x_cg(),config1_cgrange3.calc_x_cg()]
x_cg_wing_config1_range_emp = [config1_cgrange1.x_cg_wing,config1_cgrange2.x_cg_wing,config1_cgrange3.x_cg_wing]

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


