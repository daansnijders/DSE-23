# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:50:12 2019

@author: Stijn
"""
from inputs.constants import *
from inputs.concept_1 import *
import numpy as np

class empennage(object):
    def __init__(self,config,C_L_w,C_D_w,C_M_ac_w , aoa_rad):
        self.N_engine = T_req[config] * np.sin(i_e_rad)
        self.T_engine = T_req[config] * np.cos(i_e_rad)
        self.C_M_ac_w = C_M_ac_w
        self.C_N_w = C_L_w * np.cos(aoa_rad) + C_D_w * np.sin(aoa_rad)
        self.C_T_w = C_D_w * np.cos(aoa_rad) - C_L_w * np.sin(aoa_rad)
        self.config = config
        self.x_w = x_le_MAC + 0.25*MAC # 0.25?
        self.z_w = t_c * MAC /2
        
    def calc_C_N(self): 
        self.C_N = MTOW[self.config]/(0.5*rho*V_cruise**2 * S)*np.cos(theta_rad)
        return self.C_N
    
    def calc_x_h(self):
        if self.config == 1:
            A = (self.calc_C_N() - self.C_N_w - self.N_engine/(0.5*rho*V_cruise**2*S)*np.cos(i_e_rad) - self.T_engine/(0.5*rho*V_cruise**2*S)*np.sin(i_e_rad))
            print(A)
            x_h = x_cg + self.C_M_ac_w *MAC/A + self.C_N_w * (x_cg[self.config] - self.x_w)/A - self.C_T_w * (z_cg - self.z_w) + (self.T_engine * np.cos(i_e_rad) - self.N_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (z_cg - z_engine) + (self.N_engine * np.cos(i_e_rad) + self.T_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (x_cg[self.config] - x_engine)
        return x_h
    
    

#
#    
#    
#    0 = C_m_ac_w + C_N_w * (x_cg - x_w)/MAC - C_T_w * (z_cg - z_w)/MAC + calc_C_N_h() * S_h/S * (V_h/V)**2 * (x_cg - x_h)/MAC + 
#    C_N_h_times_S_h_S
#        
        
e = empennage(1,1.2,0.1,0,np.deg2rad(0.0))
print(e.calc_x_h())

"""
Values needed
M = 0
M_ac_w = sybren
M_ac_h = 0
*M_ac_c
C_N_w = sybren
C_N_h = sybren
*N_c 
T_w = tbc wrt N_w
T_h = tbc wrt N_h
*T_c
x_cg = get_cg
z_cg = get_cg
*x_w
z_w = around 0 + 0.5 * t
x_h = 0.9 of l_f?
z_h = 0.8 * d_f_outer
*x_c
*z_c
N_e = np.sin(i_e) * thrust
T_e = np.cos(i_e) * thrust
i_e = 0
x_e = get_cg
z_e = get_cg
V_cruise = concept 1
V_h = ?
V_c = ?
theta/aoa = aerodynamic guys -> alpha cruise
"""
