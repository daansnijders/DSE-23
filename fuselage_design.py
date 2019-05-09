# -*- coding: utf-8 -*-
"""
Created on Thu May  9 10:58:36 2019

@author: Sybren
"""

""" Fuselage Design """
from math import *


# [N_pax, Range [km], l_cabin, V_os, V_cc (Available), V_carry_on, V_check_in, V_cargo (V_available - V_luggage), l_f]
conc1_conf1 = [90, 4000, 19.44] 
conc1_conf2 = [120, 2000, 25.92]
conc1_conf3 = [120, 4000, 25.92]
conc2_conf1 = [90, 4000, 25.92]
conc2_conf2 = [120, 2000, 25.92]
conc2_conf3 = [120, 4000, 25.92]
conc3_conf1 = [90, 4000, 19.44]
conc3_conf2 = [120, 2000, 25.92]
conc3_conf3 = [120, 4000, 25.92]

conc_conf = [[conc1_conf1, conc1_conf2, conc1_conf3], [conc2_conf1, conc2_conf2, conc2_conf3],[conc3_conf1,conc3_conf2,conc3_conf3]]

inch_to_m = 0.0254

N_sa = 5
N_aisle = 1
aisle_width = 0.51
s_clearance = 0.02
seat_pitch = 32*inch_to_m      #[m]
seat_width = 20*inch_to_m      #[m]
armrest = 2*inch_to_m          #[m]


for i in range(len(conc_conf)):
    for j in range(len(conc_conf[i])):
        V_os = 2*0.2*conc_conf[i][j][2]*0.74
        conc_conf[i][j].append(V_os)

conc_conf[1][1][3] = conc_conf[0][1][3]
d_f_inner = N_sa*seat_width + (N_sa+N_aisle+1)*armrest + N_aisle*aisle_width + 2*s_clearance
d_f_outer = 1.045*d_f_inner + 0.084

R = d_f_outer/2
h_max = 0.9
h_floor = 0.2

p = R - h_max - h_floor
phi = 2 * acos(1-p/R)
A_cc = 0.5 * R**2 * (phi - sin(phi))
W_carry_on = 6.1                #kg
W_check_in  = 16.7              #kg
W_pax = 83.8 - W_carry_on       #kg
rho_lugg = 170                  #kg/m3
rho_cargo = 160                 #kg/m3

for i in range(len(conc_conf)):
    for j in range(len(conc_conf[i])):
        V_cc = conc_conf[i][j][2]*0.45*A_cc
        Wtot_carry_on = conc_conf[i][j][0] * W_carry_on
        Wtot_check_in = conc_conf[i][j][0] * W_check_in
        V_carry_on = Wtot_carry_on / rho_lugg
        V_check_in = Wtot_check_in / rho_lugg
        V_cargo = V_cc - (V_carry_on + V_check_in - conc_conf[i][j][3])
        M_cargo = V_cargo * rho_cargo
        M_payload = conc_conf[i][j][0] * (W_carry_on + W_check_in + W_pax) + M_cargo
        conc_conf[i][j].append(V_cc)
        conc_conf[i][j].append(V_carry_on)
        conc_conf[i][j].append(V_check_in)
        conc_conf[i][j].append(V_cargo)
        conc_conf[i][j].append(M_cargo)
        conc_conf[i][j].append(M_payload)
        print (M_payload, "kg")
        
l_nose = 1.85 * d_f_outer       #m
l_tailcone = 3.5*d_f_outer      #m
l_tail = 1.6 * d_f_outer        #m
l_cockpit = 4                   #m

for i in range(len(conc_conf)):
    for j in range(len(conc_conf[i])):
        l_f = l_cockpit + l_tail + conc_conf[i][j][2] 
        conc_conf[i][j].append(l_f)
        print(l_f)




