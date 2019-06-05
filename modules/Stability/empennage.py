# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:50:12 2019

@author: Stijn
"""
from inputs.constants import *
from inputs.concept_1 import *
from modules.Stability.Stability_runfile import x_cg_min
import numpy as np

class empennage(object):
    def __init__(self,config,C_L_w,C_D_w,C_M_ac_w , aoa_rad,x_cg):
        self.N_engine = T_req[config] * np.sin(i_e_rad)
        self.T_engine = T_req[config] * np.cos(i_e_rad)
        self.C_M_ac_w = C_M_ac_w
        self.C_N_w = C_L_w * np.cos(aoa_rad) + C_D_w * np.sin(aoa_rad)
        self.C_T_w = C_D_w * np.cos(aoa_rad) - C_L_w * np.sin(aoa_rad)
        self.config = config
        self.x_w = x_le_MAC + 0.25*MAC # 0.25?
        self.z_w = t_c * MAC /2
        self.x_cg = x_cg
        
    def calc_C_N(self): 
        self.C_N = MTOW[self.config]/(0.5*rho*V_cruise**2 * S)*np.cos(theta_rad)
        return self.C_N
    
    def calc_x_h(self):
        if self.config == 1:
            A = (self.calc_C_N() - self.C_N_w - self.N_engine/(0.5*rho*V_cruise**2*S)*np.cos(i_e_rad) - self.T_engine/(0.5*rho*V_cruise**2*S)*np.sin(i_e_rad))
            
            x_h = self.x_cg + self.C_M_ac_w *MAC/A + self.C_N_w * (self.x_cg - self.x_w)/A - self.C_T_w * (z_cg - self.z_w) + (self.T_engine * np.cos(i_e_rad) - self.N_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (z_cg - z_engine) + (self.N_engine * np.cos(i_e_rad) + self.T_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (self.x_cg - x_engine)
        return x_h
    
    
e = empennage(1,1.0,0.02,0,np.deg2rad(1.0),x_cg_min)
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
we need cl/alpha in order to determine the lift produced.
"""

class empennage2:
    def __init__(self, config, x_ac, CL_a_h, CL_a_ah, de_da, S_h, l_h, S, c, V_h, V, x_lemac):   
        self.config = config
        self.x_ac=x_ac #from nose in [m]
        self.CL_a_h = CL_a_h
        self.CL_a_ah = CL_a_ah
        self.de_da = de_da
        self.S_h = S_h
        self.l_h = l_h
        self.S = S
        self.c = c
        self.V_h = V_h
        self.V = V
        self.x_lemac = x_lemac

        self.hortail_vol = self.S_h * self.l_h / (self.S * self.c)

    
    def calc_xnp(self):
        self.x_np = (self.x_ac + self.CL_a_h / self.CL_a_ah * (1-self.de_da) * self.hortail_vol * (self.V_h / self.V)**2)*self.c
        return self.x_np
    
    def calc_xcg(self):
        self.x_cg = self.calc_xnp() - 0.05*self.c
        return self.x_cg
    
    def plot_stability(self):
        a = 1/(self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h*(self.V_h/self.V)**2)
        b = (self.x_ac - 0.05*self.c) / (self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h*(self.V_h/self.V)**2)
        
        l = np.arange(self.x_lemac, (self.x_lemac+self.c+0.01), 0.01)
        self.Sh_S = []
        for i in range (len(l)):
            Sh_S.append(a*l[i]+b)
        
        return self.Sh_S

e2 = empennage2(1, 11.78, 3.82, 4.90, 0.3835, 21.72, 16., 93.5, 3.8, 1., 1., 11.78)

r = e2.calc_xnp()    
    
