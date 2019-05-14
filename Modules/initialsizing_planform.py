# -*- coding: utf-8 -*-
"""
Created on Mon May 13 17:29:43 2019

@author: Lisa
"""
import numpy as np


def get_S(MTOW,W_S):
    return [MTOW[i] * 9.80665 /W_S[i] for i in range(3)]

def get_b(A,S):
    return [np.sqrt(A[i]*S[i]) for i in range(3)]

def get_lambda_4_rad(M_cr, M_x):
    return [np.arccos((0.75 * M_x)/(M_cr + 0.03)) for i in range(3)]

def get_taper_ratio(lambda_4_rad):
    return [0.2*(2-lambda_4_rad[i]) for i in range(3)]

def get_lambda_2_rad(lambda_4_rad,A,taper_ratio):
    return [np.arctan(np.tan(lambda_4_rad[i])-1/A[i]*(1-taper_ratio[i])/(1+taper_ratio[i])) for i in range(3)]

def get_Cr(S,taper_ratio,b):
    return [2*S[i]/((1+taper_ratio[i])*b[i]) for i in range(3)]

def get_CL(MTOW,rho,V,S):
    return [MTOW[i]*9.80665/(0.5*rho*V**2*S[i]) for i in range(3)]

def get_t_c(lambda_2_rad,M_x, M_cr,CL):
    return [(np.cos(lambda_2_rad[i])**3*(M_x-(M_cr + 0.03)*np.cos(lambda_2_rad[i]))-0.115*CL[i]**1.5)/(np.cos(lambda_2_rad[i])**2) for i in range(3)]

def get_MAC(Cr, taper_ratio):
    return [Cr[i] * 2/3 * (1+taper_ratio[i] + taper_ratio[i]**2) / (1 + taper_ratio[i]) for i in range(3)]

def get_Ct(Cr, taper_ratio):
    return [Cr[i] * taper_ratio[i] for i in range(3)]

def get_y_MAC(b, Cr, MAC, Ct):
    return [b[i]/2 * (Cr[i] - MAC[i]) / (Cr[i] - Ct[i]) for i in range(3)]

def get_dihedral_rad(lambda_4_rad):
    return [3-lambda_4_rad[i]/10 + 2 for i in range(3)]











def planformsizing(H,M_cr,S,MTOW,A):

    M_dd=M_cr+0.03
    if M_cr>0.7:
        sweep_c_4=arccos(0.75*M_x/M_dd)
    else:
        sweep_c_4=1
        
    b=sqrt(S*A)
    taper=0.2*(2-sweep_c_4*pi/180)
    
    c_r=2*S/((1+taper)*b)
    c_t=c_r*taper
    
    dihedral=deg2rad(3-rad2deg(sweep_c_4)/10) #add 2 for low wing , -2 for high wing
    
    p=p0*(1-0.0065*H/(288.15))**(9.81/(287*0.0065))
    q=0.5*1.4*p*M_cr**2
    
    Cl=MTOW/(q*S)
    MAC=c_r*2/3*(1+taper+taper**2)/(1+taper)
    Y_mac=b/2*(c_r-MAC)/(c_r-c_t)
    
    sweep_c_2=arctan(tan(sweep_c_4)-(1-taper/(a*(1+taper))))
    
    t_c=(cos(sweep_c_2)**3*(M_x-M_dd*cos(sweep_c_2))-0.115*Cl**1.5)/cos(sweep_c_2)
    if t_c>0.18:
        t_c=0.18
        
    return b,taper,c_r,c_t,dihedral,Cl, MAC, Y_mac, sweep_c_4, sweep_c_2, t_c

    
