# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""





from inputs.concept_1 import *
from inputs.constants import *

from modules.loaddiagram_detailed import *
from modules.Stability.Check_ground import *


config1     = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, OEW[0], MTOW[0], x_cg[0], y_cg[0], z_cg, MAC, S, b, A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing[0], 1)     
config2     = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, OEW[1], MTOW[1],  x_cg[1], y_cg[1], z_cg, MAC, S, b, A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing[1], 2)     
config3     = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, OEW[2], MTOW[2],  x_cg[2], y_cg[2], z_cg, MAC, S, b, A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, x_cg_wing[2], 3)     

cg1_pass = [0, 0, 0]
cg2_pass = [0, 0, 0]
weight_pass = [0, 0, 0]
cg1_pass[0], cg2_pass[0], weight_pass[0] = config1.loading_diagrams_pass()
cg1_pass[1], cg2_pass[1], weight_pass[1] = config2.loading_diagrams_pass()
cg1_pass[2], cg2_pass[2], weight_pass[2] = config3.loading_diagrams_pass()

cg1_fuel = [0, 0, 0]
cg2_fuel = [0, 0, 0]
weight_fuel = [0, 0, 0]
cg1_fuel[0], cg2_fuel[0], weight_fuel[0] = config1.loading_diagrams_fuel()
cg1_fuel[1], cg2_fuel[1], weight_fuel[1] = config2.loading_diagrams_fuel()
cg1_fuel[2], cg2_fuel[2], weight_fuel[2] = config3.loading_diagrams_fuel()


config1     = Stability_check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg[0], x_mlg[0])     
config2     = Stability_check_ground(cg1_pass[1], cg2_pass[1], weight_pass[1], cg1_fuel[1], cg2_fuel[1], weight_fuel[1], x_nlg[1], x_mlg[1])     
config3     = Stability_check_ground(cg1_pass[2], cg2_pass[2], weight_pass[2], cg1_fuel[2], cg2_fuel[2], weight_fuel[2], x_nlg[2], x_mlg[2])     

frac_min, frac_max, frac = config1.check_equilibrium()
