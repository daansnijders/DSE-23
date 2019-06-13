# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:56:29 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt

from inputs.concept_1 import *
from inputs.constants import *
from modules.main_class2 import *
from modules.Stability.cg_weight_loadingdiagram import *
from modules.Stability.cg_weight_config1 import *
from modules.Stability.control_surf_func import *
from modules.Stability.empennage import *


"""NEED FROM OTHER FILES"""
x_ac     = (x_le_MAC[0]+0.25*MAC)
CL_a_h   = 3.82
CL_a_ah  = 4.90
de_da    = 0.3835
Vh_V     = 1.
Cm_ac    = -0.3
CL_ah    = 1.6
x_cg     = x_cg_max
CL_h     = -0.5838
CL_c     = 1.3
CL_a_c   = 1.0
"""====================="""

# initialize class:
empennage1 = empennage(0, x_ac, CL_a_h, CL_a_ah, de_da, S_h, l_h[0], S, c, Vh_V, x_le_MAC[0], Cm_ac, CL_ah, x_cg, CL_h, CL_c, CL_a_c, 1., 0., 0., 0.5, 0.5, 0.5, 0.5, 1., 12.)

# outputs:
x_le_MAC = empennage1.x_le_MAC                                                  # [m] x-location of MAC main wing

S_h             = empennage1.S_h                                                # [m^2] surface area of htail
A_h             = empennage1.A_h                                                # [-] aspect ratio htail
b_h             = empennage1.b_h                                                # [m] span htail
Cr_h            = empennage1.Cr_h                                               # [m] root chord htail
Ct_h            = empennage1.Ct_h                                               # [m] tip chord htail
taper_ratio_h   = empennage1.taper_ratio_h                                      # [-] taper ratio htail
lambda_h_le_rad = empennage1.lambda_h_le_rad                                    # [rad] leading edge sweep htail
lambda_h_2_rad  = empennage1.lambda_h_2_rad                                     # [rad] half chord sweep htail
lambda_h_4_rad  = empennage1.lambda_h_4_rad                                     # [rad] quarter chord sweep htail
x_h             = empennage1.x_h                                                # [m] x-location of ac of the htail?

S_v             = empennage1.S_v                                                # [m^2] surface area of vtail
A_v             = empennage1.A_v                                                # [-] aspect ratio vtail
b_v             = empennage1.b_v                                                # [m] span vtail
Cr_v            = empennage1.Cr_v                                               # [m] root chord vtail
Ct_v            = empennage1.Ct_v                                               # [m] tip chord vtail
taper_ratio_v   = empennage1.taper_ratio_v                                      # [-] taper ratio vtail
lambda_v_le_rad = empennage1.lambda_v_le_rad                                    # [rad] leading edge sweep vtail
lambda_v_2_rad  = empennage1.lambda_v_2_rad                                     # [rad] half chord sweep vtail
lambda_v_4_rad  = empennage1.lambda_v_4_rad                                     # [rad] quarter chord sweep vtail
x_v             = empennage1.x_v                                                # [m] x-location of ac of the vtail?

# control surfaces: (inputs still need to be worked on...)
c_elev = get_c_elev(Cr_h, Ct_h, b_h)                                                       # [m] chord length elevator
S_elev = get_S_elev(S_h)                                                        # [m^2] surface area elevator
b_elev = get_b_elev(S_elev,c_elev)                                              # [m] span elevator

c_rud = get_c_rud(Cr_v)                                                         # [m] chord length rudder
S_rud = get_S_rud(S_v)                                                          # [m^2] surface area rudder
b_rud = get_b_rud(S_rud,c_rud)                                                  # [m] span rudder

c_ail = get_c_ail(Ct)                                                           # [m] chord length aileron
S_ail = get_S_ail(S)                                                            # [m^2] surface area aileron
b_ail = get_b_ail(b)                                                            # [m] span aileron

c_splr = get_c_splr(MAC)                                                        # [m] chord length spoiler
b_splr = get_b_splr(b)                                                          # [m] span spoiler