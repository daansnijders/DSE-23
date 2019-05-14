# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:44:09 2019

@author: thong

Analysis of Concept 1 fuselage joint.
120pax, 4000km

Assumptions:
1) Fuselage is cylindrical
2) Fuselage structural mass and payload is distributed equally
3) Clamped at wing
4) Tail is not generating lift
5) z_CG at the center of fuselage, symmetrical in x,y,zplane
6) Aluminium 7075-T6, yield strength 503MPa
"""
import numpy as np
import matplotlib.pyplot as plt
from isa import isa
from mpl_toolkits.mplot3d import Axes3D

"""
Parameters
Previous update time: - 
Current update time: 13/5/2019, 12:11PM
-------------------------------------------------------------------------------
"""
# MTOW [kg] @ 90,4000km
MTOW1 = 49391.28866
# MTOW [kg] @ 120,4000km
MTOW2 = 68264.27835
# OEW [kg]
OEW = 38729.81
# Safety factor
safety_factor = 1.5
# Load Factor
load_factor = 2.1
# Gravitation acceleration [m/s^2]
g = 9.81

# Fuselage length [m]
l_fuselage = 29.335        
# Fuselage diameter [m]
d_fuselage = 3.685
# Fuselage structure mass [kg]
m_fuselage = 15973.84113

# Payload mass [kg]
m_payload = 12925.76943

# Wing lift [N]
L_wing = MTOW1*g
# Wing group mass [kg]
m_winggroup = 11946.24871
# Wing group location [m]
x_wing = 16.1349902
# Fuel mass [kg]
m_fuel = 16608.69892


# Tail mass [kg]
m_tail = 1638.34268
# Tail location [m]
x_tail = 26.4015


# Canard mass [kg]
m_canard = (157.4050988-113.8873926)/157.4050988*6280.313608
# Canard location [m]
x_canard = 8

"""
Adjusting m_fuselage to include main and nose landing gear
"""
delta_oew = OEW - (m_fuselage+m_winggroup+m_canard+m_tail)
m_fuselage+=delta_oew


""" Force and Moment Equilibrium """
delta_F = (m_fuselage+m_payload+m_winggroup+m_canard+m_tail+m_fuel)*g-L_wing
delta_M = -m_fuselage*l_fuselage/2*g-m_payload*l_fuselage/2*g-\
            m_winggroup*x_wing*g-m_canard*g*x_canard-m_tail*g*x_tail-m_fuel*g*x_wing+L_wing*x_wing

# Canard lift [N]
L_canard = (-delta_M-x_tail*delta_F)/(x_canard-x_tail)
# Tail lift [N]
L_tail = delta_F-L_canard



"""Analysis"""
n = 1000
step_size = l_fuselage/n

"""
Point Loads: Canard weight, Canard Lift,Fuel Weight, Wing Group weight, Wing Lift, Tail weight, Tail Lift
Distributed Loads: Fuselage, Payload
"""
P_loads = np.array([-m_canard*g, L_canard, -m_fuel*g, -m_winggroup*g, L_wing, -m_tail*g, L_tail])*load_factor
P_x = np.array([x_canard,x_canard,x_wing,x_wing,x_wing,x_tail,x_tail])
D_loads = np.array([-m_fuselage*g/n, -m_payload*g/n])*load_factor
D_x = np.linspace(0,l_fuselage,n)

"""
Internal Shear Diagram
"""

x = 0
x_lst = []
shear_lst = []
shear = 0
for i in range(n):
    x += step_size
    force = 0
    for j in range(len(P_loads)):
        if x-step_size/2<P_x[j]<=x+step_size/2:
            shear += P_loads[j]
    shear += sum(D_loads)
    x_lst.append(x)
    shear_lst.append(shear)


plt.figure("Internal Shear")
plt.plot(x_lst,shear_lst)
plt.title("Internal Shear Diagram")
plt.xlabel("x [m]")
plt.ylabel("Internal Shear [Nm]")
plt.show()

"""
Internal Moment Diagram
"""
moment_lst = []
moment = 0
for i in range(len(x_lst)):
    if i < 1:
        moment += (shear_lst[i] + shear_lst[i+1])/2*step_size
    else:
        moment += (shear_lst[i] + shear_lst[i-1])/2*step_size
    moment_lst.append(moment)


plt.figure("Internal Moment Diagram")
plt.plot(x_lst,(moment_lst))
plt.title("Internal Moment Diagram")
plt.xlabel("x [m]")
plt.ylabel("Internal Moment [Nm]")
plt.show()

""" Fuselage cross-section design"""
# Fuselage Stringer Area [mm^2]
stringer_area = 120
# Tensile Yield Strength [MPa]
yield_strength = 503

"""Calculate MMOI of fuselage"""
def mmoi(stringer_number,radius,stringer_area):   #[-,mm,mm^2]
    angle = 2*np.pi/stringer_number
    stringer_lst = []
    for i in range(stringer_number):
        y = np.sin(angle*i)*radius
        z = np.cos(angle*i)*radius
        stringer_lst.append([y,z,stringer_area])
    stringer_lst = np.array(stringer_lst)
    mean_y = sum(stringer_lst[:,2]*stringer_lst[:,1])/sum(stringer_lst[:,2])
    MMOI = sum((stringer_lst[:,1])**2*stringer_lst[:,2])
    return [MMOI,mean_y]

"""Pressure at 2438ft [Pa]"""
P1 = isa(2438/3.281)[1]
"""Pressure at 37000ft [Pa]"""
P2 = isa(37000/3.281)[1]

"""Calculate number of stringers needed"""
fuselage_radius = d_fuselage/2     #[m]
stringer_number = 2
stringer_lst = []
stress_lst = []

i = 0
while i < len(moment_lst):
    k = moment_lst[i]
    long_stress = (P1-P2)*np.pi*(fuselage_radius)**2/(stringer_number*stringer_area)
    bending_stress = (abs(k)*1000*fuselage_radius*1000/(mmoi(stringer_number,fuselage_radius*1000,stringer_area)[0]))
    normal_stress = ((bending_stress)+long_stress)*safety_factor
    
    if yield_strength<normal_stress:
        stringer_number += 1
    else:
        stringer_lst.append(stringer_number)
        stress_lst.append(normal_stress)
        i += 1
        stringer_number = 1
        

plt.figure("Number of Stringers")
plt.plot(x_lst,stringer_lst)
plt.title("Number of Stringers")
plt.xlabel("x location(m)")
plt.ylabel("Number of Stringers")
plt.show()

"""Number of bolts at the joint"""

x_joint1 = -(25.92-19.44)/2+x_canard         #[m]
x_joint2 = (25.92-19.44)/2+x_canard         #[m]
n_joint1 = stringer_lst[min(range(len(x_lst)), key=lambda i: abs(x_lst[i]-x_joint1))]
n_joint2 = stringer_lst[min(range(len(x_lst)), key=lambda i: abs(x_lst[i]-x_joint2))]

"""Bolts Details"""
l_reinforcement = 500           # [mm]
rho_material = 2810             # [kg/m^3]


""" Thickness of skin """
t_skin = (P1-P2)*safety_factor*fuselage_radius/(yield_strength*10**6)
print("Skin thickness= ",t_skin*10**3, " mm")

""" Maximum Stress """
print("Maximum Normal Stress = ", max(stress_lst), " MPa")
    
""" Shear Stress """
def shear(shear_lst, stringer_lst,fuselage_radius, stringer_area):
    radius = fuselage_radius*1000
    shear_stress = []
    for i in range(len(shear_lst)):
        stringer = []
        for j in range(stringer_lst[i]):
            angle = 2*np.pi/stringer_lst[i]
            y = np.cos(angle*j)*radius
            z = np.sin(angle*j)*radius
            Q = radius*radius*np.cos(angle*j) #<------ thickness removed
            q = - shear_lst[i] * safety_factor*Q/mmoi(stringer_lst[i],radius,stringer_area)[0]
            stringer = [shear_lst[i],x_lst[i],y,z,Q,q]
            shear_stress.append(stringer)
    shear_stress = np.array(shear_stress)
    return shear_stress

test = shear(shear_lst, stringer_lst,fuselage_radius, stringer_area)
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(test[:,1],test[:,2],test[:,3],c=abs(test[:,-1]))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
    
    