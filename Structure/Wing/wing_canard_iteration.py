# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 11:59:33 2019

@author: thong
"""

import Structure.Wing.lift_surface_function as lsf
import numpy as np
import inputs.constants as c

rho = c.rho_0

def wing_struc_analysis(CL_max,velocity_max,span,root_chord,tip_chord,surface_area):
    
    """Obtain Lift at every point of the wing"""
    lift = 0.5*rho*CL_max*surface_area*velocity_max**2
    q_lift = lift/surface_area

    N = 100
    step_size = span/(2*N)
    force_lst = []
    y_lst = []
    for i in range(N):
        yi = i*step_size
        cri = lsf.get_chord(yi,span,root_chord,tip_chord)
        cti = lsf.get_chord(yi+step_size,span,root_chord,tip_chord)
        Si = 0.5*(cri+cti)*step_size
        Li = q_lift*Si
        force_lst.append(Li)
        y_lst.append(yi)
    force_lst = np.array(force_lst)
    y_lst = np.array(y_lst)
    
    """Obtain internal bending moment at every point of the wing"""
    root_moment = sum(force_lst*y_lst)
    root_shear = -sum(force_lst)
    int_bending_moment_lst = []
    for i in range(N):
        moment_i = 0
        y = y_lst[i]
        for j in range(len(force_lst)):
            if y_lst[j]<yi:
                moment_i += force_lst[j]*(y-y_lst[j])
        moment = root_moment + moment_i + root_shear*y
        int_bending_moment_lst.append(moment)
    int_bending_moment_lst = np.array(int_bending_moment_lst)

    """ Obtain stringer number """
    l_wing_box = 0.5
    h_wing_box = (0.1102+0.1083)/2
    t_skin = 0.003
    t_spar = 0.015
    
    stringer_area = 400*10**-6
    strength = 344*10**6
    density = 2800
    geo_prop_lst = []
    
    for i in range(N):
        top_stringer = 0
        bottom_stringer = 0
        moment = int_bending_moment_lst[i]
        y = y_lst[i]
        c = lsf.get_chord(y,span,root_chord,tip_chord)
        height = h_wing_box*c
        length = l_wing_box*c
        
        moi_skin = 2*((1/12)*length*t_skin**3+length*t_skin*(height/2)**2)
        moi_spar = 2*((1/12)*height*t_spar**3)
        moi = moi_skin+moi_spar
        
        stress = moment*height/(2*moi)
        
        while stress>strength:
            top_stringer += 1
            bottom_stringer += 1
            moi_stringer = (top_stringer+bottom_stringer)*stringer_area*\
                            (height/2)**2
            moi_total = moi+moi_stringer
            stress = moment*height/(2*moi_total)
            
        enclosed_area = length*height
        material_area = 2*(length*t_skin+height*t_spar)+stringer_area*(top_stringer+\
                          bottom_stringer)
        fuel_area = enclosed_area-material_area
        geo_prop_lst.append([enclosed_area,material_area,fuel_area,top_stringer,bottom_stringer])
    geo_prop_lst = np.array(geo_prop_lst)
    
    """ Obtain mass """
    total_mass = density*sum(geo_prop_lst[:,1]*step_size)*2
    fuel_mass_available = sum(geo_prop_lst[:int(3/4*N),2]*step_size)*804*2
    
    """ Obtain total cost """
    total_cost = total_mass*4.37*2
    
    return fuel_mass_available

#print(wing_struc_analysis(1.584,240,35.44,5.698,1.763,132.211))
        
        
        
    
    
    