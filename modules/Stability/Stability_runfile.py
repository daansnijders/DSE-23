# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:43:09 2019

@author: daansnijders
"""


from modules.class2_struct_defs import *
from modules.class2_powerplant_defs import *
from modules.class2_fixedequipment_defs import *
from modules.class2 import *

from inputs.concept_1 import *
from inputs.constants import *
from inputs.performance_inputs import *

from modules.loaddiagram_detailed import *


config1     = Loading_diagram(x_cargo[0], l_f[0], l_cabin[0], seat_pitch, N_pax[0], N_sa, OEW[0], MTOW[0], x_cg[0], y_cg[0], z_cg[0], MAC[0], S[0], b[0], A, Xfirst, M_payload[0], M_cargo_available[0], M_fuel[0], M_pax, M_carry_on, x_cg_wing[0], 1)     
config2     = Loading_diagram(x_cargo[1], l_f[1], l_cabin[1], seat_pitch, N_pax[1], N_sa, OEW[1], MTOW[1],  x_cg[1], y_cg[1], z_cg[1], MAC[1], S[1], b[1], A, Xfirst, M_payload[1], M_cargo_available[1], M_fuel[1], M_pax, M_carry_on, x_cg_wing[1], 2)     
config3     = Loading_diagram(x_cargo[2], l_f[2], l_cabin[2], seat_pitch, N_pax[2], N_sa, OEW[2], MTOW[2],  x_cg[2], y_cg[2], z_cg[2], MAC[2], S[2], b[2], A, Xfirst, M_payload[2], M_cargo_available[2], M_fuel[2], M_pax, M_carry_on, x_cg_wing[2], 3)     

cg_min = [0, 0, 0]
cg_max = [0, 0, 0]
weight_max = [0, 0, 0]
cg_min[0], cg_max[0], weight_max = config1.loading_diagrams()
cg_min[1], cg_max[1], weight_max = config2.loading_diagrams()
cg_min[2], cg_max[2], weight_max = config3.loading_diagrams()
