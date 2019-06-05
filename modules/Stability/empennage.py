# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:50:12 2019

@author: Stijn
"""
from inputs.constants import *
import numpy as np

class empennage(object):
    def __init__(self,C_N_w,C_N_c,V,V_h,V_c,S,S_h,S_c,N_e,T_e,i_e_rad,theta_rad,W,config):
        self.C_N_w = C_N_w
        self.C_N_c = C_N_c
        self.V = V
        self.V_h = V_h
        self.V_c = V_c
        self.S = S
        self.S_h = S_h
        self.S_c = S_c
        self.N_e = N_e
        self.T_e = T_e
        self.i_e_rad = i_e_rad
        self.theta_rad = theta_rad
        self.W = W
        self.config = config
        
    def calc_C_N(self):
        self.C_N = self.W/(0.5*rho*self.V**2 * self.S)*np.cos(self.theta_rad)
        return self.C_N
    
    def calc_C_N_h(self):
        if self.config == 1:
            self.C_N_h = (self.calc_C_N() - self.C_N_w - self.N_e/(0.5*rho*self.V**2*self.S)*np.cos(self.i_e_rad) - self.T_e/(0.5*rho*self.V**2*self.S)*np.sin(self.i_e_rad))/((self.V_h/self.V)**2*self.S_h/self.S)

        return self.C_N_h
    
    C_N_w = C_L * np.cos(aoa_rad) + C_D * np.sin(aoa_rad)
    C_T_w = C_D * np.cos(aoa_rad) - C_L * np.sin(aoa_rad)

    
    
    0 = C_m_ac_w + C_N_w * (x_cg - x_w)/MAC - C_T_w * (z_cg - z_w)/MAC + calc_C_N_h() * S_h/S * (V_h/V)**2 * (x_cg - x_h)/MAC + (T_e * np.cos(i_e_rad)\
                            - N_e * np.sin(i_e_rad)) / (0.5*rho*V**2*S*MAC) * (z_cg - z_e) + (N_e * np.cos(i_e_rad) + T_e * np.sin(i_e_rad)) / (0.5*rho*V**2*S*MAC) * (x_cg - x_e)
    C_N_h_times_S_h_S
        
        
e = empennage(1,3,4,5,6,7,8,9,10,11,12,13,14,1)

"""
Values needed
M = 0
M_ac_w = aerodynamic guys
M_ac_h = 0
*M_ac_c
N_w = tbc
N_h = tbc
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
theta/aoa = aerodynamic guys -> alpha cruise
"""
