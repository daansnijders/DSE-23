# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:43:50 2019

@author: thong
"""

import numpy as np
import lift_surface_function as lsf

def get_fuel_force(cross_section_lst,N,wing_span,load_factor,fuel_mass,LE_sweep,front_spar,rear_spar,Cr,Ct):
    """
    Calculate amount of fuel it can carry 
    """
    step_size = (wing_span/2)/N
    tot_fuel_mass = 0
    for i in range(len(cross_section_lst)):
        if cross_section_lst[i,0]<=0.75*wing_span/2 and cross_section_lst[i,0]>=0.1*wing_span/2:
            tot_fuel_mass += cross_section_lst[i,3]*step_size*804*2
    
    print("Total available fuel mass [kg] = ",tot_fuel_mass)
    
    """
    Calculate fuel weight at every point
    [xi,yi,zi,Fx_i,Fy_i,Fz_i]
    """
    fuel_force_lst = []
    fuel = fuel_mass
    for i in range(N):
        yi = i*step_size
        zi = 0
        C_i = lsf.get_chord(yi,wing_span,Cr,Ct)        
        xi = yi*np.tan(LE_sweep/180*np.pi)+cross_section_lst[i,4]
        fuel_i = (cross_section_lst[i,3]*(wing_span/2)/N)*804
        fuel -= fuel_i
        if fuel>0:
            fuel_force_lst.append([xi,yi,zi,0,0,fuel_i*-9.81*load_factor])
        else:
            fuel_force_lst.append([xi,yi,zi,0,0,0])

    return(tot_fuel_mass,fuel_force_lst)