# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:44:09 2019

@author: moosh

Analysis of Concept 1 fuselage joint.
120pax, 4000km

Assumptions:
1) Fuselage is cylindrical
2) Fuselage structural mass and payload is distributed equally
3) Point load, unclamped, fuel and wing weight conincide with wing lift, tail lift
    conincide with tail weight, canard lift conincide with canard weight
4) Stringers carry longitudinal stress
5) z_CG at the center of fuselage, symmetrical in x,y,zplane
6) Aluminium 7075-T6, yield strength 503MPa
7) Torsion neglected
8) Effect of buckling neglected
9) Calculation all done to meet limit load according the tensile yield strength

-------------------------------------------------------------------------------
"""


def fueslage_stress_max(l_fuselage, x_cg_hwing, x_cg_vwing, x_cg_ngear, x_cg_mgear, x_cg_wing_group, Xstart, Xlast, X_wingbox_start,\
X_wingbox_end, M_payload, M_fuselage, M_fittings, M_horizontal_tail, M_vertical_tail, M_landinggear_nose, M_landinggear_main, M_wing_group,\
Lift_mainwing, Lift_tail):
    import numpy as np
    import matplotlib.pyplot as plt
    from isa import isa
    
    n= 1000
    step_size = l_fuselage/n
    stringer_area = 0.0004
    stringer_no = [0] *n
    for i in range(n):
        stringer_no[i]=60
        
    radius = [0] *n
    for i in range(n):
        radius[i]=3.685
    
    MOI = [0] *n
    for i in range(n):
        MOI = stringer_no[i]*(radius[i])**2 *stringer_area * 0.5
    
    payload_w = [0] *n
    for i in range(int((Xfirst/lfuselage)*n),int((Xlast/lfuselage)*n)):
        payload_w[i]=-M_payload*9.80665/(int((Xlast/lfuselage)*n)-int((Xfirst/lfuselage)*n)) 
        
    fuselage_w = [0] *n
    for i in range(n):
        fuselage_w[i]=-(M_fuselage+M_fittings)*9.80665/n
        
    component_w = [0] *n 
    component_w[int((x_cg_hwing/lfuselage)*n)]=-M_horizontaltail*9.80665
    component_w[int((x_cg_vwing/lfuselage)*n)]=-M_verticaltail*9.80665
    component_w[int((x_cg_ngear/lfuselage)*n)]=-M_landinggear_nose*9.80665
    component_w[int((x_cg_mgear/lfuselage)*n)]=-M_landinggear_main*9.80665
    component_w[int((x_cg_wing_group/lfuselage)*n)]=-M_wing_group*9.80665
       
    
    
    lift_forces = [0] *n 
    for i in range(int((X_wingbox_start/lfuselage)*n),int((X_wingbox_end/lfuselage)*n)):
        lift_forces[i]=Lift_mainwing/(int((X_wingbox_end/lfuselage)*n)-int((X_wingbox_start/lfuselage)*n))
        
    lift_forces[int((x_cg_hwing/lfuselage)*n)]=Lift_tail
    
    forces_sum =  [0] *n 
    for i in range(n):
        foces_sum[i]= payload_w[i]+fuselage_w[i]+component_w[i]+lift_forces[i]
    
    V =  [0] *n 
    for i in range(n):
        if i == 0:
            V[i]=forces_sum[i]
        else:
            V[i]=V[i-1]+forces_sum[i]
    
    M =  [0] *n 
    for i in range(n):
        for j in range(i):
            M[i]+=forces_sum[j]*step_size*j
    
    p_diff = isa(2438/3.281)[1] - isa(37000/3.281)[1]
    stress_pressure_long =  [0] *n
    for i in range(n):
        stress_pressure_long[i]=(p_diff*np.pi*(radius[i])**2)/(stringer_area*stringer_no[i])
        
    stress_bending_max = [0] *n
    for i in range(n):
        stress_bending_max[i] = M[i]*radius[i]/MOI[i]
        
    stress_long_max = [0] *n
    for i in range(n):
        stress_long_max[i] = abs(stress_bending_max[i])+stress_pressure_long[i]
        
        return max(stress_bending_max)