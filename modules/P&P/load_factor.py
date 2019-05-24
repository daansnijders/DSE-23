# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:16:35 2019

@author: Stijn, Rik
"""
import numpy as np
from inputs.constants import *
from inputs.concept_1 import * 

# inputs:
C_L_max = 1.5                                                                   # [-]
C_D_C_L_max = 0.05                                                              # [-]
MTOW = 65000                                                                    # [kg]
S = max(S)

def calc_V_S(MTOW,C_L_max,C_D_C_L_max,S,rho):
    W = MTOW *9.80665                                                           # [N]
    C_N_max = np.sqrt(C_L_max**2 + C_D_C_L_max**2)                              # [-]
    V_S = np.sqrt(2*W/(S * rho * C_N_max))                                      # [m/s]
    return V_S

calc_V_S(MTOW,C_L_max,C_D_C_L_max,S,rho)

def calc_n_lim_pos(MTOW):
    W = MTOW * 2.2046226218488
    n_lim_pos = 2.1 + 24000 / (W + 10000)
    return n_lim_pos

calc_n_lim_pos(MTOW)
