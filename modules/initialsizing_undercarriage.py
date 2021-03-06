# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:30:05 2019

@author: Stijn
"""
import numpy as np                                                  

def get_P_mw(MTOW,N_mw,weight_distribution):
    return [(1-weight_distribution)*MTOW[i]*9.80665/N_mw for i in range(3)]                 

def get_P_nw(MTOW,N_nw,weight_distribution):
    return [weight_distribution*MTOW[i]*9.80665/N_nw for i in range(3)]

def get_wheel_sizing_mw(tire_pressure,P_mw):
    return [(tire_pressure/100,P_mw[i]/9.80665) for i in range(3)]

def get_wheel_sizing_nw(tire_pressure,P_nw):
    return [(tire_pressure/100,P_nw[i]/9.80665) for i in range(3)]

def tangent(x,y,angle_rad):
    return np.tan(angle_rad),y- np.tan(angle_rad) * x

def get_x_mlg(z_cg,theta_rad, beta_rad, x_cg, stroke, l_f):
    beta_rad_correct = beta_rad- np.pi/2
    
    scrape = [tangent(l_f[i],z_cg,theta_rad) for i in range(3)]
    tip_over = [tangent(x_cg[i],z_cg,beta_rad_correct) for i in range(3)]
    
    x_mlg = [(tip_over[i][1] - scrape[i][1] + stroke)/(scrape[i][0] - (tip_over[i][0])) for i in range(3)] 
    return x_mlg


def get_z_mlg(x_mlg,beta_rad,x_cg, z_cg):
    beta_rad_correct = beta_rad - np.pi/2

    tip_over = [tangent(x_cg[i],z_cg,beta_rad_correct) for i in range(3)]
    
    
    z_mlg = [tip_over[i][0] * x_mlg[i] + tip_over[i][1] for i in range(3)]
    return z_mlg

def get_l_mw(x_mlg,x_cg):
    return [x_mlg[i] - x_cg[i] for i in range(3)]

def get_z_mlg_new(l_mw,beta_rad,z_cg):
    beta_rad_correct = beta_rad - np.pi/2
    z_mlg=[np.tan(beta_rad_correct)*(l_mw[i])+z_cg for i in range(3)]
    return z_mlg


def get_new_beta_config23(z_mlg,z_cg,l_mw):
    z_tot=z_cg-z_mlg
    beta_new_rad=np.pi/2-np.arctan(z_tot/l_mw)
    return beta_new_rad

def check_new_scrape_angle(l_f,x_mlg,d_f_outer,z_mlg):
    theta_new_rad=[np.arctan((d_f_outer-z_mlg)/(l_f[i]-x_mlg[i])) for i in range(3)]
    return theta_new_rad
    
def get_l_nw(l_w,P_mw,N_mw,P_nw,N_nw):
    return [(l_w[i]*P_mw[i]*N_mw)/(P_nw[i]*N_nw) for i in range(3)]

def get_x_nlg(x_cg,l_n):
    return [x_cg[i] - l_n[i] for i in range(3)]

def get_y_mlg(b,dihedral_rad,psi_rad,phi_rad,z_cg,z_mlg,l_n,l_w,y_eng,z_eng,d_eng):
    z_t = b/(2) * np.tan(dihedral_rad)
    
    y_mlg1 = [(l_n[i] + l_w[i])/(np.sqrt((l_n[i]**2 \
              *np.tan(psi_rad)**2)/(z_cg - z_mlg)**2 - 1)) for i in range(3)]
    y_mlg2 = [y_eng - (z_eng - z_mlg)/np.tan(phi_rad) for i in range(3)]                           
    y_mlg3 = [b/2 - (z_t-z_mlg)/np.tan(phi_rad) for i in range(3)] 
    
    y_mlg = [max([y_mlg1[i],y_mlg2[i],y_mlg3[i]]) for i in range(3)]                                      
    return y_mlg

# Get diameter of landing gear strut
# Input: force on landing gear [N], length of strut [m]
# Output: Diameter of landing gear [m]
def get_d_lg(force,length):
    d = 0.01                       #[m]
    #Material: Low alloy steel, 300M, quenched and tempered
    strength = 1590 * 10**6         #[Pa]
    comp_stress = force/(np.pi*(d/2)**2)
    bend_stress = (0.8*force*length*d/2)/(np.pi*(d**4)/64)
    stress = (comp_stress+bend_stress)*1.5
    while stress>strength:
        comp_stress = force/(np.pi*(d/2)**2)
        bend_stress = (0.8*force*length*d/2)/(np.pi*(d**4)/64)
        stress = (comp_stress+bend_stress)*1.5
        d += 0.01
    return d

