# -*- coding: utf-8 -*-
"""
Created on Thu May  9 10:58:36 2019

@author: Sybren
"""

""" Fuselage Design """
from math import *

#Initial data for the different concepts and configurations
# [N_pax, Range [km], l_cabin, V_os, V_cc (Available), V_carry_on, V_check_in, V_cargo (V_available - V_luggage), l_f]
conc1_conf1 = [94, 4000, 19.44] 
conc1_conf2 = [125, 2000, 25.92]
conc1_conf3 = [125, 4000, 25.92]
conc2_conf1 = [94, 4000, 19.44]
conc2_conf2 = [125, 2000, 25.92]
conc2_conf3 = [125, 4000, 25.92]
conc3_conf1 = [94, 4000, 19.44]
conc3_conf2 = [125, 2000, 25.92]
conc3_conf3 = [125, 4000, 25.92]

#generate one list
conc_conf = [[conc1_conf1, conc1_conf2, conc1_conf3], [conc2_conf1, conc2_conf2, conc2_conf3],[conc3_conf1,conc3_conf2,conc3_conf3]]

#Go from inches to meters
inch_to_m = 0.0254

#data on the interior
N_sa = 5                        # Number of seats abreast
N_aisle = 1                     # Number of aisles
aisle_width = 0.51              #[m]
s_clearance = 0.02              #[m] Clearance between seat and fuselage
seat_pitch = 32*inch_to_m       #[m]
seat_width = 20*inch_to_m       #[m]
armrest = 2*inch_to_m           #[m] Armrest width


for i in range(len(conc_conf)):
    for j in range(len(conc_conf[i])):
        V_os = 2*0.2*conc_conf[i][j][2]*0.74        #[m3] Overhead storage volume
        conc_conf[i][j].append(V_os)
        
conc_conf[1][0][2] = 25.92      #[m] The length of the fuselage does not change, but module does not contain overhead space

d_f_inner = N_sa*seat_width + (N_sa+N_aisle+1)*armrest + N_aisle*aisle_width + 2*s_clearance    #[m] Inner diameter
d_f_outer = 1.045*d_f_inner + 0.084                                                             #[m] Outer diameter

R = d_f_outer/2     #[m] Radius of the fuselage
h_max = 0.9         #[m] Distance between the center of the fuselage and the floor
h_floor = 0.2       #[m] Floor heigth

p = R - h_max - h_floor     #[m] Distance between the lower point of the inner fuselage and the floor
phi = 2 * acos(1-p/R)       #[rad] Angle between the two connection points of the fuselage and floor
A_cc = 0.5 * R**2 * (phi - sin(phi)) #[Cargo hold area]
W_carry_on = 6.1                #[kg] Average carry-on luggage weight
W_check_in  = 16.7              #[kg] Average check-in luggage weight
W_pax = 83.8 - W_carry_on       #[kg] Average passenger weight
rho_lugg = 170                  #[kg/m3] Average luggage density
rho_cargo = 160                 #[kg/m3] Average cargo density

for i in range(len(conc_conf)):
    for j in range(len(conc_conf[i])):
        V_cc = conc_conf[i][j][2]*0.45*A_cc                 #[m3] Total available cargo hold volume
        Wtot_carry_on = conc_conf[i][j][0] * W_carry_on     #[kg] Total carry-on weight
        Wtot_check_in = conc_conf[i][j][0] * W_check_in     #[kg] Total check-in weight
        print (Wtot_carry_on + Wtot_check_in)
        V_carry_on = Wtot_carry_on / rho_lugg               #[m3] Total carry-on volume needed
        V_check_in = Wtot_check_in / rho_lugg               #[m3] Total check-in volume needed
        V_cargo = V_cc - (V_carry_on + V_check_in - conc_conf[i][j][3]) #[m3] Total available cargo volume
        M_cargo = V_cargo * rho_cargo                       #[kg] Available cargo weight
        M_payload = conc_conf[i][j][0] * (W_carry_on + W_check_in + W_pax) + M_cargo #[kg] Total payload weight
        conc_conf[i][j].append(V_cc)
        conc_conf[i][j].append(V_carry_on)
        conc_conf[i][j].append(V_check_in)
        conc_conf[i][j].append(V_cargo)
        conc_conf[i][j].append(M_cargo)
        conc_conf[i][j].append(M_payload)
        
l_cockpit = 4                   #[m] Cockpit length
l_nose = 1.85 * d_f_outer       #[m] Nose cone length
l_tailcone = 3.5*d_f_outer      #[m] Tail cone length
l_tail = 1.6 * d_f_outer        #[m] Tail length

for i in range(len(conc_conf)):
    for j in range(len(conc_conf[i])):
        l_f = l_cockpit + l_tail + conc_conf[i][j][2]       #[m]
        conc_conf[i][j].append(l_f)



