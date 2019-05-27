# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:25:39 2019

@author: thong

Analysis of the mass contribution of the rail system that connects the fuselage
and the wing.

Assumptions:
1) Single rail system
2) Mass of cabin distributed equally across the rail which has a length of short cabin length
3) No mass at structure above the module
    
Descriptions:

"""

import numpy as np
from isa import isa
import matplotlib.pyplot as plt
"""
RAIL ANALYSIS PARAMETERS
Previous update time: - 
Current update time: 15/5/2019, 12:11PM
-------------------------------------------------------------------------------
"""
# Largest Fuselage mass [kg]
fuselage_mass = 16528.09546
# Largest Fuselage length [m]
fuselage_length = 35.815
# Largest Payload mass [kg]
payload_mass = 13339.56529
# Largest Cabin Length [m]
cabin_length = 25.92

# Fuselage Diameter [m]
fuselage_diameter = 3.685

# Shortest Cabin Length [m]
short_cabin_length = 19.44

# tail_length [m]
tail_length = 5.895

# Cabin height [m]
cabin_height = 2.5

# Largest module mass [kg]
module_mass = (cabin_length+tail_length)/fuselage_length*fuselage_mass

# Gravitational Acceleration [m/s^2]
g = 9.81

# Safety factor 
safety_factor = 1.5

# Load factor
load_factor = 2.5

# Density of material [kg/m^3]
rho = 2810
# Tensile Yield Strength [MPa]
strength = 159

# Bank Angle
bank_angle = 30/180*np.pi

# Total Force [N]
total_force = (module_mass+payload_mass)*load_factor*safety_factor*g

""" 
Rail Design Parameters
-------------------------------------------------------------------------------
Wing Rail
"""
b1 = 45              #[mm]
t1 = 10              #[mm]

"""
Fuselage Rail
"""
b2 = 20             #[mm]
b3 = 45             #[mm]
t2 = 10             #[mm]


"""
FUSELAGE DESIGN PARAMETERS
Previous update time: - 
Current update time: 13/5/2019, 12:11PM
-------------------------------------------------------------------------------
"""

# MTOW [kg] @ 120,4000km
MTOW = 70632.89
# OEW [kg]
OEW = 40108.34


# Payload mass [kg]
m_payload = 13339.56786

# Wing lift [N]
L_wing = MTOW*g
# Wing group mass [kg]
m_winggroup = 12360.75575
# Wing group location [m]
x_wing = 14.18011683               #16.1349902
# Fuel mass [kg]
m_fuel = 17184.98214


# Tail mass [kg]
m_tail = 1695.18936
# Tail location [m]
x_tail = 32.2335            #26.4015

# Nose and Main Landing Gear Location [m]
x_nlg = 4
x_mlg = 19

"""Pressure at 2438ft [Pa]"""
P1 = isa(2438/3.281)[1]
"""Pressure at 37000ft [Pa]"""
P2 = isa(37000/3.281)[1]

""" Fuselage cross-section design"""
# Fuselage Stringer Area [mm^2]
stringer_area = 400

"""
Analysis (Step Size)
-------------------------------------------------------------------------------
"""

n = 1000
step_size = fuselage_length/n


""" 
Wing Rail Analysis ****All stress calculated in MPa**** [VERIFY!!!!!!!!]
-------------------------------------------------------------------------------
Calculate the stress of the rail of the wing side based on tensile stress
"""
tensile_stress_w = total_force/(b1*short_cabin_length*1000)

"""
Calculate the bending stress of rail of the aircraft at banking
"""
bending_stress_w = (cabin_height/2*1000+t1+t2)*total_force*np.sin(bank_angle)*b1/2 \
                    /(1/12*short_cabin_length*1000*b1**3)   #(mm+mm+mm)*N*-*mm/(mm*mm^3)
             
"""
Difference in strength and induced stress
"""
delta_normal_stress = strength-(bending_stress_w+tensile_stress_w)

"""
Bending and shear stress at rail teeth
"""
bending_stress_teeth_w = (total_force/2*b2/2)*t1/2/(1/12*short_cabin_length*1000*t1**3)
shear_stress_teeth_w = (total_force)/(short_cabin_length*1000*t1)

"""
Fuselage Rail Analysis ****All stress calculated in MPa****
-------------------------------------------------------------------------------
Calculate the tensile stress of the section between the fuselage and teeth
"""
tensile_stress_t_f = total_force/(2*b3*short_cabin_length*1000)

"""
Bending stress at the connection between the fuselage and teeth due to banking of aircraft
"""
bending_stress_t_f = (cabin_height/2*1000+t1/2)*total_force*np.sin(bank_angle)*b3/2 \
                    /(1/12*short_cabin_length*1000*b3**3)   #(mm+mm+mm)*N*-*mm/(mm*mm^3)

"""
Bending and shear stress at rail teeth
"""
bending_stress_teeth_f = (total_force/2*b2/2)*t2/2/(1/12*short_cabin_length*1000*t2**3)
shear_stress_teeth_f = (total_force)/(short_cabin_length*1000*t2)

"""
Mass of rail
"""
rail_volume = (b1+b2*2+b3*2)*(t1+t2)*short_cabin_length
rail_mass = (b1+b2*2+b3*2)*(t1+t2)*10**-6*rho*short_cabin_length

print("---------------------Rail Analysis-----------------------------")
print("Material strength = ", strength, " MPa")
print(" ")
print("Tensile Stress at wing rail = ", tensile_stress_w," MPa")
print("Bending Stress at wing rail due to banking = ", bending_stress_w," MPa")
print("Bending Stress at wing rail teeth = ", bending_stress_teeth_w," MPa")
print("Shear Stress at wing rail teeth = ", shear_stress_teeth_w," MPa")
print()
print("Tensile Stress at fuselage rail teeth connection = ", tensile_stress_t_f," MPa")
print("Bending Stress at fuselage rail teeth connection due to banking = ", bending_stress_t_f," MPa")
print("Bending Stress at fuselage rail teeth = ", bending_stress_teeth_f," MPa")
print("Shear Stress at fuselage rail teeth = ", shear_stress_teeth_f," MPa")

print("")
print("Rail mass", rail_mass, " kg")
print("Rail width = ", b1+b2*4, "mm" )
print("Rail thickness = ", t1+t2, "mm" )

"""
Fuselage analysis
-------------------------------------------------------------------------------
"""
print(" ")
"""
-------------------------------------------------------------------------------
Adjusting m_fuselage to include main and nose landing gear
"""
delta_oew = OEW - (fuselage_mass+m_winggroup+m_tail)
print("Delta OEW = ", delta_oew, " kg")
fuselage_mass += delta_oew/2
m_nlg = delta_oew/8
m_mlg = 3/8*delta_oew

""" 
-------------------------------------------------------------------------------
Force and Moment Equilibrium 
in flight 
-------------------------------------------------------------------------------
"""
delta1_F = (fuselage_mass+m_payload+m_winggroup+m_tail+m_fuel+m_mlg+m_nlg)*g
delta1_M = -fuselage_mass*fuselage_length/2*g-m_payload*fuselage_length/2*g-\
            m_winggroup*x_wing*g-m_tail*g*x_tail-m_fuel*g*x_wing-m_mlg*g*x_mlg\
            -m_nlg*g*x_nlg
            
print ("Delta force" , delta1_F)
print ("Delta Moment" , delta1_M)

# Wing lift [N]
L_wing = (-delta1_M-x_tail*delta1_F)/(x_wing-x_tail)
# Tail lift [N]
L_tail = delta1_F-L_wing
print(MTOW*g-L_wing)

EQ1_F = L_wing+L_tail-(fuselage_mass+m_payload+m_winggroup+m_tail+m_fuel+m_mlg+m_nlg)*g
EQ1_M = L_wing*x_wing+L_tail*x_tail\
        -(fuselage_mass*fuselage_length/2+m_payload*fuselage_length/2+\
          m_winggroup*x_wing+m_tail*x_tail+m_fuel*x_wing+m_nlg*x_nlg+m_mlg*x_mlg)*g
print("Summation Force = ", EQ1_F)
print("Summation Moment = ", EQ1_M)

"""
In flight
Point Loads: Fuel Weight, Wing Group weight, Wing Lift, Tail weight, Tail Lift
Distributed Loads: Fuselage without Module , Module, Payload
"""
P1_loads = np.array([-m_nlg*g,-m_mlg*g,-m_fuel*g, -m_winggroup*g, L_wing, -m_tail*g, L_tail])*load_factor
P1_x = np.array([x_nlg,x_mlg, x_wing,x_wing,x_wing,x_tail,x_tail])
D1_loads = np.array([-(fuselage_mass)*g/fuselage_length, -m_payload*g/fuselage_length])*load_factor
D1_x = np.linspace(0,fuselage_length,n)

"""
Internal Shear Diagram
"""

x = 0
x_lst = []
shear1_lst = []
shear1 = 0
for i in range(n):
    x += step_size
    for j in range(len(P1_loads)):
        if x-step_size/2<P1_x[j]<=x+step_size/2:
            shear1 += P1_loads[j]
    shear1 += sum(D1_loads)*step_size
    x_lst.append(x)
    shear1_lst.append(shear1)


plt.figure("Internal Shear")
plt.plot(x_lst,shear1_lst)
plt.title("Internal Shear Diagram")
plt.xlabel("x [m]")
plt.ylabel("Internal Shear [Nm]")
plt.show()

"""
Internal Moment Diagram
Macaulay's Step Function
"""
moment1_lst = []
moment1 = 0
for i in range(n):
    x_loc = i*step_size
    moment1 = sum(D1_loads)*0.5*x_loc**2
    number = 0
    for j in range(len(P1_x)):
        if x_loc>P1_x[j]:
            moment1 +=  (x_loc-P1_x[j])*P1_loads[j]
            number += 1

    moment1_lst.append(moment1)
print("End internal moment = ", moment1_lst[-1])


plt.figure("Internal Moment Diagram")
plt.plot(x_lst,(moment1_lst))
plt.title("Internal Moment Diagram")
plt.xlabel("x [m]")
plt.ylabel("Internal Moment [Nm]")
plt.show()

""" 
-------------------------------------------------------------------------------
Force and Moment Equilibrium 
on ground
"""
delta2_F = ((fuselage_mass-module_mass)+m_winggroup+m_tail+m_fuel+m_mlg+m_nlg)*g
delta2_M = -(fuselage_mass-module_mass)*fuselage_length/2*g-\
            m_winggroup*x_wing*g-m_tail*g*x_tail-m_fuel*g*x_wing-m_mlg*g*x_mlg\
            -m_nlg*g*x_nlg
            
print ("Delta force2" , delta2_F)
print ("Delta Moment2" , delta2_M)

# Force from Nose Landing Gear [N]
F_nlg = (-delta2_M-x_mlg*delta2_F)/(x_nlg-x_mlg)
# Force from Main Landing Gear [N]
F_mlg = delta2_F-F_nlg


EQ2_F = F_nlg+F_mlg-((fuselage_mass-module_mass)+m_winggroup+m_tail+m_fuel+m_mlg+m_nlg)*g
EQ2_M = F_nlg*x_nlg+F_mlg*x_mlg\
        -((fuselage_mass-module_mass)*fuselage_length/2+\
          m_winggroup*x_wing+m_tail*x_tail+m_fuel*x_wing+m_nlg*x_nlg+m_mlg*x_mlg)*g
print("Summation Force2 = ", EQ2_F)
print("Summation Moment2 = ", EQ2_M)

"""
In ground
Point Loads: Fuel Weight, Wing Group weight, Wing Lift, Tail weight, Tail Lift
Distributed Loads: Fuselage without Module , Module, Payload
"""
P2_loads = np.array([-m_nlg*g,-m_mlg*g,-m_fuel*g, -m_winggroup*g, F_nlg, -m_tail*g, F_mlg])
P2_x = np.array([x_nlg,x_mlg, x_wing,x_wing,x_nlg,x_tail,x_mlg])
D2_loads = np.array([-(fuselage_mass-module_mass)*g/fuselage_length])
D2_x = np.linspace(0,fuselage_length,n)

"""
Internal Shear Diagram
"""

x = 0
x_lst = []
shear2_lst = []
shear2 = 0
for i in range(n):
    x += step_size
    for j in range(len(P2_loads)):
        if x-step_size/2<P2_x[j]<=x+step_size/2:
            shear2 += P2_loads[j]
    shear2 += sum(D2_loads)*step_size
    x_lst.append(x)
    shear2_lst.append(shear2)


plt.figure("Internal Shear")
plt.plot(x_lst,shear2_lst)
plt.title("Internal Shear Diagram")
plt.xlabel("x [m]")
plt.ylabel("Internal Shear [Nm]")
plt.show()

"""
Internal Moment Diagram
Macaulay's Step Function
"""
moment2_lst = []
moment2 = 0
for i in range(n):
    x_loc = i*step_size
    moment2 = sum(D2_loads)*0.5*x_loc**2
    number = 0
    for j in range(len(P2_x)):
        if x_loc>P2_x[j]:
            moment2 +=  (x_loc-P2_x[j])*P2_loads[j]
            number += 1

    moment2_lst.append(moment2)
print("End internal moment = ", moment2_lst[-1])


plt.figure("Internal Moment Diagram")
plt.plot(x_lst,(moment2_lst))
plt.title("Internal Moment Diagram")
plt.xlabel("x [m]")
plt.ylabel("Internal Moment [Nm]")
plt.show()

"""
Moment of Inertia of the entire Fuselage
-------------------------------------------------------------------------------
"""
def mmoi1(t,fuselage_diameter,cabin_height):                 #[mm,m,m]
    area_fuselage = fuselage_diameter*1000*np.pi*t              #[mm^2]
    mmoi_fuselage = np.pi*(fuselage_diameter*1000/2)**3*t       #[mm^4]
    
    alpha = np.arccos((fuselage_diameter/2-(fuselage_diameter-cabin_height))/(fuselage_diameter/2))
    beam_width = 2*1000*(fuselage_diameter/2)*np.sin(alpha) #[mm]
    z_beam = (fuselage_diameter/2-(fuselage_diameter-cabin_height))*1000 #[mm]
    area_beam = beam_width*t #[mm^2]
    mmoi_beam = (1/12)*beam_width*t**3         #[mm^4]
    z_neutral = (area_beam*z_beam)/(area_beam+area_fuselage)     #[mm]
    area = area_beam+area_fuselage
    mmoi = mmoi_beam + mmoi_fuselage + area_fuselage*z_neutral**2 + area_beam*(z_beam-z_neutral)**2
    max_z = fuselage_diameter/2*1000-z_neutral #[mm]
    return [max_z,z_neutral,area,mmoi]

"""
Moment of Inertia of one section of fuselage
-------------------------------------------------------------------------------
"""
def mmoi2(t,fuselage_diameter,cabin_height):            #[mm,m,m]
    alpha = np.arccos((fuselage_diameter/2-(fuselage_diameter-cabin_height))/(fuselage_diameter/2))
    beam_width = 2*1000*(fuselage_diameter/2)*np.sin(alpha) #[mm]
    z_beam = (fuselage_diameter/2-(fuselage_diameter-cabin_height))*1000 #[mm]
    area_beam = beam_width*t #[mm^2]
    mmoi_beam = (1/12)*beam_width*t**3         #[mm^4]
    """Discretise the arc and obtain moment of inertia"""
    arc_length = fuselage_diameter*1000/2*2*alpha       #[mm]
    arc_area = arc_length*t                             #[mm^2]
    dis_n = 1000
    angle = np.linspace(-alpha,alpha,dis_n)
    arc_z = fuselage_diameter/2*1000*np.cos(angle)  #[mm]
    disc_area = arc_area/dis_n #[mm^2]
    """Area and mmoi of entire structure except cabin"""
    z_neutral = (area_beam*z_beam + sum(disc_area*arc_z))/(arc_area+area_beam)   #[mm]
    area = arc_area+area_beam       #[mm^2]
    mmoi = sum(disc_area*(arc_z-z_neutral)**2) + mmoi_beam+area_beam*(z_beam-z_neutral)**2 #[mm^2]
    max_z = (fuselage_diameter/2)*1000 - z_neutral #[mm]
    return [max_z,z_neutral,area,mmoi]



"""
Analysis of required thickness
"""
with_module = []
with_module_t = []
without_module = []
without_module_t = []

"""
Thickness at the flight without module
"""
t = 0.1
j = 0
while j < len(moment1_lst):
    I = mmoi2(t,fuselage_diameter,cabin_height)
    moment = moment1_lst[j]*1000        #[Nmm]
    stress = abs(moment*I[0])*safety_factor/I[3]
    if strength<stress:
        t += 0.1
    else:
        without_module.append([I[0],I[1],I[2],I[3],t,stress])
        without_module_t.append(t)
        t=0.1
        j+=1

without_module = np.array(without_module)
"""
Thickness at the flight with module
"""
T = 0.1
j = 0
while j < len(moment1_lst):
    I = mmoi1(t,fuselage_diameter,cabin_height)
    moment = moment1_lst[j]*1000        #[Nmm]
    stress = abs(moment*I[0])*safety_factor/I[3]
    if strength<stress:
        t += 0.1
    else:
        with_module.append([I[0],I[1],I[2],I[3],t,stress])
        with_module_t.append(t)
        t=0.1
        j+=1

with_module = np.array(with_module)

"""
Calculate reinforcement mass needed
"""

volume_diff = 0
volume_withmodule = 0
check1 = []
check2 = []
for i in range(len(with_module_t)):
    mmoi_1 = mmoi2(without_module_t[i],fuselage_diameter,cabin_height) #without module with thickness of without module
    mmoi_2 = mmoi2(with_module_t[i],fuselage_diameter,cabin_height) # without module with thickness of with module
    mmoi_3 = mmoi1(with_module_t[i],fuselage_diameter,cabin_height) # Entire structure with thickness of with module
    check1.append(mmoi_1)
    check2.append(mmoi_2)
    area_difference = mmoi_1[2]-mmoi_2[2]
    volume_withmodule += mmoi_3[2]*step_size*1000        #[mm^3]
    volume_diff += area_difference*step_size*1000        #[mm^3]

check1 = np.array(check1)
check2 = np.array(check2)
mass_reinforcement = rho*volume_diff*10**-9
percentage = (volume_diff)/volume_withmodule*100

print(" ")
print("[ANS]   Percentage differece in fuselage mass= ",percentage, " %")
print("[ANS]   Fuselage mass= ",mass_reinforcement, " kg")

#print(mmoi2(1,fuselage_diameter,cabin_height)[-1]-mmoi1(1,fuselage_diameter,cabin_height)[-1])
#
#print(mmoi2(1,fuselage_diameter,cabin_height))
#print(mmoi1(1,fuselage_diameter,cabin_height))
plt.figure("Fuselage Thickness")
plt.plot(x_lst,without_module_t,label="Remaining Fuselage")
plt.plot(x_lst,with_module_t,label="Fuselage with Cabin")
plt.title("Fuselage Skin Thickness")
plt.legend()
plt.xlabel("x [m]")
plt.ylabel("Skin Thickness [mm]")
plt.show()