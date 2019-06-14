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
import matplotlib.pyplot as plt
import numpy as np
import geo_prop_function as gpf
import fuel_force_function as fff

import sys
sys.path.insert(0,'/Users/thong/Documents/TU Delft DSE/DSE-23')
import inputs.concept_1 as c1

"""
Parameters
------------------------------------------------------------------------------
"""
# Front Spar location [-]
spar_front = gpf.spar_front

# Rear Spar location [-]
spar_rear = gpf.spar_rear

# Maximum Take-off Mass [kg]
MTOM = c1.MTOW[0]
print("MTOM [kg] = ",MTOM)

# Wing span [m]
wing_span = c1.b
print("Wing span [m] = ",wing_span)

# Root chord [m]
chord_root = c1.Cr
print("Root chord [m] = ",chord_root)

# Tip Chord [m]
chord_tip = c1.Ct
print("Tip chord [m] = ",chord_tip)
# Wing surface area [m^2]
wing_area = c1.S
print("Wing area [m^2] = ",wing_area)

# Discretisation number [-]
N = 200

# Lift Coefficient [-]
coeff_lift = 1.584
# Drag Coefficient [-]
coeff_drag = 0.00668

# Load_factor [-]
load_factor = 3.1

# Safety_factor [-]
safety_factor = 1.5

# Leading edge sweep [deg]
sweep_LE = c1.lambda_le_rad*180/np.pi

# Engine property list[x_loc(m), y_loc(m), z_loc(m), thrust [N], y_force, weight[N]
engine_thrust = 108540    #[N]
engine_mass = 2177      #[kg]
print("Y engine [m] = ",c1.y_engine) 
engine = [0,c1.y_engine,0,-engine_thrust,0,engine_mass*-9.81*load_factor]

# main landing gear property list[x_loc(m), y_loc(m), z_loc(m), drag [N], y_force, weight[N]]
mlg_drag = 200       #[N]
mlg_mass = 2500     #[kg]
x_mlg = c1.x_mlg[0]-c1.x_le_w[0]
y_mlg = c1.y_mlg[0]
z_mlg = c1.z_mlg
print("X mlg [m] = ",x_mlg)
print("Y mlg [m] = ",y_mlg)
print("Z mlg [m] = ",z_mlg)
mlg = [x_mlg,y_mlg,-z_mlg,mlg_drag,0,mlg_mass*-9.81*load_factor]

# Fuel Mass [kg]
fuel_mass = c1.M_fuel[0]
print("fuel mass [kg] = ",fuel_mass)

"""
Form list of cross-section
"""
cross_section_lst = []
for i in range(N):
    y_loc = i*wing_span/(2*N)
    chord = lsf.get_chord(y_loc,wing_span,chord_root,chord_tip)
    cross_section = gpf.get_cross_sec_prop(chord,y_loc,wing_span)
    cross_section_lst.append(cross_section)
cross_section_lst = np.array(cross_section_lst)

"Calculate fuel force"
fuel = fff.get_fuel_force(cross_section_lst,N,wing_span,load_factor,\
                                   fuel_mass,sweep_LE,spar_front,spar_rear,\
                                   chord_root,chord_tip)
tot_fuel_mass = fuel[0]
fuel_force_lst = fuel[1]

check = sum(np.array(fuel_force_lst)[:,-1])/load_factor/9.81
print(check)

"""
Calculate structural weight force
"""
struc_force_lst = gpf.get_struc_force(cross_section_lst,load_factor,sweep_LE,wing_span,\
                                      chord_root,chord_tip)
 
"""
Calculate total cost and mass of wing box
"""
mass_cost = gpf.get_wingbox_mass_cost(cross_section_lst,wing_span)
print("Wing box mass [kg] = ", mass_cost[0])
print("Wing box cost [USD] = ", mass_cost[1])
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
force_lst.extend(fuel_force_lst)
force_lst.extend(struc_force_lst)
force_lst = np.array(force_lst)
force_lst=force_lst[np.argsort(force_lst[:,1])]

#plt.figure()
#plt.scatter(force_lst[:,1],force_lst[:,0])

"""
Form list of internal loads
"""
internal_loads_lst = lsf.get_internal_force_distri(wing_span, chord_root,\
                                                   chord_tip,N,force_lst,\
                                                   sweep_LE,spar_front,spar_rear)
internal_loads_lst = np.array(internal_loads_lst)

"""
ANSWER ANALYSIS 1
"""
#plt.figure("Shear x")
#plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,3])
#plt.figure("Shear z")
#plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,4])
#plt.figure("Moment x")
#plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,5])
#plt.figure("Moment z")
#plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,6])
#plt.figure("Torque")
#plt.plot(internal_loads_lst[:,0],internal_loads_lst[:,7])
#plt.show()

"""
Calculate normal stress due to bending
sigma = -My/I
"""
sigmax_upper = -safety_factor*internal_loads_lst[:,5]*cross_section_lst[:,7]/cross_section_lst[:,8]
sigmax_lower = -safety_factor*internal_loads_lst[:,5]*cross_section_lst[:,6]/cross_section_lst[:,8]

"""
Calculate normal stress due to engine
sigma = -Mx/I
"""
sigmaz_upper = -safety_factor*internal_loads_lst[:,6]*cross_section_lst[:,14]/cross_section_lst[:,12]
sigmaz_lower = -safety_factor*internal_loads_lst[:,6]*cross_section_lst[:,15]/cross_section_lst[:,12]

"""
Calculate critical shear web buckling stress [Pa]
[ribs_pitch,spar_height,a/b,ks_value,tau_cr]
"""
tau_cr_web_lst = []
for i in range(len(cross_section_lst)):
    tau_cr = gpf.get_web_buckling(cross_section_lst[i,0],cross_section_lst[i,9],wing_span)
    tau_cr_web_lst.append(tau_cr)
tau_cr_web_lst = np.array(tau_cr_web_lst)

"""
Calculate critical skin buckle at upper skin
"""
sigma_cr_skin_lst = []
for i in range(len(cross_section_lst)):
    sigma_cr = gpf.get_skin_buckling(cross_section_lst[i,0],cross_section_lst[i,11],wing_span)
    sigma_cr_skin_lst.append(sigma_cr)
sigma_cr_skin_lst = np.array(sigma_cr_skin_lst)

"""
Calculate shear flow at web due to lift bending
q = VQ/I
"""
q_web_lst = np.abs(internal_loads_lst[:,4]*cross_section_lst[:,10]/(cross_section_lst[:,8]))

"""
Calculate shear flow at skin due to engine
q = VQ/I
"""
q_skin_lst = np.abs(internal_loads_lst[:,3]*cross_section_lst[:,13]/(cross_section_lst[:,12]))

"""
Calculate shear due to torque
q = Torque/2A
"""
q_torque_lst = np.abs(internal_loads_lst[:,7]/(2*cross_section_lst[:,1]))

"""Calculate shear at the web"""
shear_web_lst = safety_factor*(q_torque_lst+q_web_lst)/gpf.t_spar

"""Calculate shear at the skin"""
shear_skin_lst = safety_factor*(q_torque_lst+q_skin_lst)/gpf.t_skin

"""Combine stresses"""
stress = abs(sigmax_upper)+abs(sigmaz_upper)

plt.figure("Skin Stress")
plt.plot(cross_section_lst[:,0],stress)
plt.plot(cross_section_lst[:,0],sigma_cr_skin_lst[:,4:])
plt.show

#plt.figure("Stress")
#plt.plot(cross_section_lst[:,0],sigmax_upper)
#plt.plot(cross_section_lst[:,0],sigmaz_upper)
#plt.show()