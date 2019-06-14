# -*- coding: utf-8 -*-
"""
Created on Mon May 13 17:29:43 2019

@author: Lisa
"""
import numpy as np

def get_S(MTOW,W_S):
    return [MTOW[i] * 9.80665 /W_S[i] for i in range(3)]

def get_b(A,S):
    assert np.sqrt(A*S) < 36
    return np.sqrt(A*S)

def get_lambda_4_rad(M_cr, M_x):
    return np.arccos((0.75 * M_x)/(M_cr + 0.03))

def get_lambda_4_rad_from_lambda_le(lambda_le_rad,Cr,b,taper_ratio):
    lambda_4_rad= [np.arctan(np.tan(lambda_le_rad)+Cr[i]/(2*b[i])*(taper_ratio-1)) for i in range(3)]
    return lambda_4_rad

def get_taper_ratio(lambda_4_rad):
    return 0.2*(2-lambda_4_rad)

def get_lambda_2_rad(lambda_4_rad,A,taper_ratio):
    return np.arctan(np.tan(lambda_4_rad)-1/A*(1-taper_ratio)/(1+taper_ratio))

def get_lambda_2_rad_canard(lambda_4_rad,A,taper_ratio):
    return np.arctan(np.tan(lambda_4_rad)-1/A*(1-taper_ratio)/(1+taper_ratio))

def get_Cr(S,taper_ratio,b):
    return 2*S/((1+taper_ratio)*b)

def get_Cr_canard(S,taper_ratio,b):
    return 2*S/((1+taper_ratio)*b)

def get_lambda_le_rad(lambda_4_rad, Cr, b, taper_ratio):
    return np.arctan(np.tan(lambda_4_rad)-(Cr/(2*b))*(taper_ratio-1)) 

def get_CL(MTOW,rho,V,S):
    return [MTOW[i]*9.80665/(0.5*rho*V**2*S) for i in range(3)]
def get_CL_canard(MTOW,rho,V,S):
    return [MTOW[i]*9.80665/(0.5*rho*V**2*S[i]) for i in range(1,3)]

def get_t_c(lambda_2_rad,M_x, M_cr,CL):
    return (np.cos(lambda_2_rad)**3*(M_x-(M_cr + 0.03)*np.cos(lambda_2_rad))-0.115*CL**1.5)/(np.cos(lambda_2_rad)**2)

def get_t_c_canard(lambda_2_rad,M_x, M_cr,CL):
    return [(np.cos(lambda_2_rad)**3*(M_x-(M_cr + 0.03)*np.cos(lambda_2_rad))-0.115*CL[i]**1.5)/(np.cos(lambda_2_rad)**2) for i in range(1,3)]


def get_MAC(Cr, taper_ratio):
    return Cr * 2/3 * (1+taper_ratio + taper_ratio**2) / (1 + taper_ratio)


def get_MAC_canard(Cr, taper_ratio):
    return [Cr[i] * 2/3 * (1+taper_ratio + taper_ratio**2) / (1 + taper_ratio) for i in range(1,3)]

def get_Ct(Cr, taper_ratio):
    return Cr * taper_ratio 

def get_Ct_canard(Cr, taper_ratio):
    return Cr * taper_ratio

def get_y_MAC(b, Cr, MAC, Ct):
    return b/2 * (Cr - MAC) / (Cr - Ct)

def get_y_MAC_canard(b, Cr, MAC, Ct):
    return [b[i]/2 * (Cr[i] - MAC[i]) / (Cr[i] - Ct[i]) for i in range(1,3)]

def get_dihedral_rad(lambda_4_rad):
    return np.deg2rad(3-np.rad2deg(lambda_4_rad)/10 + 2)

def get_le_wing(y_MAC,x_le_MAC, lambda_2_rad, MAC, Cr):
    return (x_le_MAC + 0.5 * MAC) - np.tan(lambda_2_rad) * y_MAC - 0.5*Cr
    
