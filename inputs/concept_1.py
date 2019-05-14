# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""
from modules.initialsizing_weights import *

T_W    = [0.295,0.295,0.295]                                                    # [-]
W_S    = [4253, 4253, 4253]                                                     # [N/m^2]
M_ff   = [0.7567, 0.8274, 0.7567]                                               # [kg]
OEW = [27745.73, 27745.73, 38729.81]                                           # [-]

MTOW = get_MTOW(OEW)
M_fuel = get_M_fuel(MTOW,M_ff)
T_req = get_T_req(T_W, MTOW)
M_payload = get_M_payload(MTOW,OEW,M_fuel)