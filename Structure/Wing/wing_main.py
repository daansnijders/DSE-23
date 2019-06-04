# -*- coding: utf-8 -*-
"""
Created on Tue May 28 14:05:25 2019

@author: thong
Wing Analysis 1 Main File

Description:
All units in N,kg,m,pa

Assumptions:
1) Point load from engine and main landing gear
------------------------------------------------------------------------------
"""
import lift_surface_function as lsf
import geo_prop_function as gpf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
Parameters
------------------------------------------------------------------------------
"""
# Maximum Take-off Mass [kg]
MTOM = 58722.6

# Wing span [m]
wing_span = 38.593274

# Root chord [m]
chord_root = 5.359098

# Tip Chord [m]
chord_tip = 1.657860

# Wing surface area [m^2]
wing_area = 135.403711

# Discretisation number [-]
N = 1000

# Lift Coefficient [-]
coeff_lift = 1.584
# Drag Coefficient [-]
coeff_drag = 0.00668

# Load_factor [-]
load_factor = 3.1

# Safety_factor [-]
safety_factor = 1.5

# Leading edge sweep [deg]
sweep_LE = 15

# Engine property list[x_loc(m), y_loc(m), z_loc(m), thrust [N], y_force, weight[N]
engine_thrust = 1000    #[N]
engine_mass = 3800      #[kg] 
engine = [0,5,-1,-engine_thrust,0,engine_mass*-9.81*load_factor]

# main landing gear property list[x_loc(m), y_loc(m), z_loc(m), drag [N], y_force, weight[N]]
mlg_drag = 200       #[N]
mlg_mass = 2500     #[kg]
mlg = [6,3,-1,mlg_drag,0,mlg_mass*-9.81*load_factor]

"""
Form list of external forces
[x,y,z,Fx,Fy,Fz]
-------------------------------------------------------------------------------
"""

force_lst = []
aero_force_lst = lsf.get_aeroforce_distri(MTOM, load_factor, wing_area, \
                                          wing_span,sweep_LE,chord_root, \
                                          chord_tip,N,coeff_lift,coeff_drag)
force_lst.extend(aero_force_lst)
force_lst.extend([mlg])
force_lst.extend([engine])
force_lst = np.array(force_lst)

"""
Form list of internal loads
"""
internal_loads_lst = lsf.get_internal_force_distri(wing_span, chord_root, chord_tip,N,force_lst,sweep_LE)
internal_loads_lst = np.array(internal_loads_lst)

"""
ANSWER ANALYSIS 1
"""
plt.figure("Shear x")
plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,3])
plt.figure("Shear z")
plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,4])
plt.figure("Moment x")
plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,5])
plt.figure("Moment y")
plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,6])
plt.figure("Torque")
plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,7])
plt.show()

"""
Form list of cross-section
"""
cross_section_lst = []
for i in range(len(internal_loads_lst)):
    cross_section = gpf.get_cross_sec_prop(internal_loads_lst[i,1],internal_loads_lst[i,0],wing_span)
    cross_section_lst.append(cross_section)
cross_section_lst = np.array(cross_section_lst)

"""
Calculate amount of fuel it can carry
"""
fuel_volume = 0
for i in range(len(cross_section_lst)):
    if internal_loads_lst[i,0]<0.6*wing_span/2:
        fuel_volume += cross_section_lst[i,2]*(wing_span/2)/N
fuel_mass = 2*fuel_volume*804