# -*- coding: utf-8 -*-
"""
Created on Tue May 28 14:10:50 2019

@author: thong
"""

import numpy as np

"""
Get chord at the y_location on span
"""
def get_chord(point,span,Cr,Ct):
    chord = Cr-(Cr-Ct)/(span/2)*point
    return chord
"""
Get aerodynamic force at every point
"""
def get_aeroforce_distri(mass,load_factor,S,b,LE_sweep,Cr,Ct,number,CL,CD,CM):
    print("------------------------------------------")
    print("Run get_aeroforce_distri")
    aero_force_lst = []
    step_size = b/(2*number)
    lift_q = mass*9.81*load_factor/S
    drag_q = mass*9.81*CD/(CL*S)
    lift = 0
    area = 0
    for i  in range(number):
        yi = i*step_size
        zi = 0
        Ct_i = get_chord(yi+step_size,b,Cr,Ct)
        Cr_i = get_chord(yi,b,Cr,Ct)        
        xi = yi*np.tan(LE_sweep/180*np.pi)+0.25*Cr_i
        Si = 0.5*(Cr_i+Ct_i)*step_size
        #Fy_i = 0
        Fx_i = drag_q*Si
        Fz_i = lift_q*Si
        T_i = lift_q*CM/CL*Cr_i
        aero_force_lst.append([xi,yi,zi,Fx_i,T_i,Fz_i])
        lift += Fz_i
        area += Si
    delta_lift = 2*lift-mass*9.81*load_factor
    delta_area = area*2-S
    print("-----Checking-----")
    print("Difference in lift [N]= ",delta_lift)
    print("Difference in area [m^2]= ",delta_area)
    return aero_force_lst
"""
Get internal force at every point
[y, chord, surface area, shear_x, shear_z, moment_x, moment_z, torque]
"""
def get_internal_force_distri(span,Cr,Ct,N,force_lst,LE_sweep,spar_front,spar_rear):
    print("------------------------------------------")
    print("Run get_internal_force_distri")
    step_size = span/(2*N)
    root_shear_z = -sum(force_lst[:,-1])
    root_shear_x = -sum(force_lst[:,3])
    root_moment_z = sum(force_lst[:,3]*force_lst[:,1])
    root_moment_x = sum(force_lst[:,-1]*force_lst[:,1])
    torque_lst = force_lst[:,5]*(force_lst[:,0]-(force_lst[:,1]*np.tan(LE_sweep/180*np.pi)+\
                      (spar_front+spar_rear)/2*get_chord(force_lst[:,1],span,Cr,Ct)))+ \
                          force_lst[:,3]*-force_lst[:,2] - force_lst[:,4]
    root_torque = sum(torque_lst)
    internal_loads = []
    print("-----Checking-----")
    print("Root Shear x[N] =  ", root_shear_x)
    print("Root Shear z[N] =  ", root_shear_z)
    print("Root Moment x[Nm] =  ", root_moment_x)
    print("Root Moment z[Nm] =  ", root_moment_z)
    print("Root Torque [Nm] =  ", root_torque)
    for i in range(N):
        yi = i*step_size
        chord = get_chord(yi,span,Cr,Ct)
        Ct_i = get_chord(yi+step_size,span,Cr,Ct)
        Cr_i = chord        
        Si = 0.5*(Cr_i+Ct_i)*step_size
        force_lst=force_lst[np.argsort(force_lst[:, 1])]
        shear_z_i = 0
        shear_x_i = 0
        moment_x_i = 0
        moment_z_i = 0
        torque_i = 0
        for j in range(len(force_lst)):
            if yi>=force_lst[j,1]:
                shear_z_i += force_lst[j,-1]
                moment_x_i += force_lst[j,-1]*(yi-force_lst[j,1])
                torque_i += torque_lst[j]
                shear_x_i += force_lst[j,3]
                moment_z_i += force_lst[j,3]*(yi-force_lst[j,1])
        torque = root_torque - torque_i
        moment_x = root_moment_x+moment_x_i+root_shear_z*yi
        shear_z = root_shear_z + shear_z_i
        shear_x = root_shear_x + shear_x_i
        moment_z = root_moment_z + moment_z_i + root_shear_x*yi
        internal_loads.append([yi,chord,Si,shear_x,shear_z,moment_x,moment_z,torque])
    print("Tip Shear x[N] =  ", internal_loads[-1][3])
    print("Tip Shear z[N] =  ", internal_loads[-1][4])
    print("Tip Moment x[Nm] =  ", internal_loads[-1][5])    
    print("Tip Moment z[Nm] =  ", internal_loads[-1][6])    
    print("Tip Torque [Nm] =  ", internal_loads[-1][7])
    return internal_loads


