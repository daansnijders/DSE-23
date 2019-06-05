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
    
    def C_N_h(self):
        if self.config == 1:
            self.C_N_h = (self.calc_C_N() - self.C_N_w - self.N_e/(0.5*rho*self.V**2*self.S)*np.cos(self.i_e_rad) - self.T_e/(0.5*rho*self.V**2*self.S)*np.sin(self.i_e_rad))/((self.V_h/self.V)**2*self.S_h/self.S)
        else:
            pass
        return self.C_N_h
        
        
e = empennage(1,3,4,5,6,7,8,9,10,11,12,13,14,1)

"""
Values needed
M = 0
M_ac_w
M_ac_h
M_ac_c
N_w
N_h
N_c
T_w
T_h
T_c
x_cg
z_cg
x_w
z_w
x_h
z_h
x_c
z_c
N_e
T_e
i_e
x_e
z_e
"""
