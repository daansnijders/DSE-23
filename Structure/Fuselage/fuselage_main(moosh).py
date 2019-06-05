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
#import sys
#sys.path.insert(0,'/Users/thong/Documents/TU Delft DSE/DSE-23')
#import inputs.concept_1 as c1
#import modules.CG.CG_func as cgfunc
#import modules.CG.class2_CG as cl2cg
#l_fuselage = c1.l_f[2]
#x_cg_hwing = cl2cg.get_cg_hwing

#def fuselage_stress_max(l_fuselage, x_cg_hwing, x_cg_vwing, x_cg_ngear, x_cg_mgear, x_cg_wing_group, Xfirst, Xlast, X_wingbox_start,\
#X_wingbox_end, M_payload, M_fuselage, M_fittings, M_horizontaltail, M_verticaltail, M_landinggear_nose, M_landinggear_main, M_wing_group,\
#Lift_mainwing, Lift_tail):

import numpy as np
import matplotlib.pyplot as plt
from isa import isa
import sys
sys.path.insert(0,'/Users/thong/Documents/TU Delft DSE/DSE-23')
import inputs.concept_1 as c1
import modules.CG.CG_func as cgfunc
import modules.CG.class2_CG as cl2cg

from modules.class2 import *
from main import *



x_cg_wing_group = 16
X_wingbox_start = 15
X_wingbox_end = 17.5
M_payload = 12500
M_fuselage = 9000
M_fittings = 9000
M_verticaltail = 1000
M_landinggear_nose = 200
M_landinggear_main = 1500
M_wing_group = 12000
Lift_mainwing = 500000
Lift_tail = -10000

l_fuselage = config1_class2.l_f
x_cg_hwing = c1.x_cg_tail
x_cg_vwing = c1.x_cg_tail
x_cg_ngear = c1.x_nlg
x_cg_mgear = c1.x_mlg
Xfirst = c1.l_cockpit
Xlast = Xfirst+c1.l_cabin
M_horizontal_tail = config1_class2.M_horizontaltail


n= 1000
step_size = float(l_fuselage)/n
x = np.linspace(0,l_fuselage,n)

stringer_area = 0.0004
stringer_no = [0] *n
for i in range(n):
    stringer_no[i]=60

radius = [0] *n
for i in range(n):
    radius[i]=3.685

MOI = [0] *n
for i in range(n):
    MOI[i] = stringer_no[i]*(radius[i])**2 *stringer_area * 0.5

payload_w = [0] *n
for i in range(int((Xfirst/l_fuselage)*n),int((Xlast/l_fuselage)*n)):
    payload_w[i]=-M_payload*9.80665/(int((Xlast/l_fuselage)*n)-int((Xfirst/l_fuselage)*n)) 
    
fuselage_w = [0] *n
for i in range(n):
    fuselage_w[i]=-(M_fuselage+M_fittings)*9.80665/n
    
component_w = [0] *n 
component_w[int((x_cg_hwing/l_fuselage)*n)]=-M_horizontaltail*9.80665
component_w[int((x_cg_vwing/l_fuselage)*n)]=-M_verticaltail*9.80665
component_w[int((x_cg_ngear/l_fuselage)*n)]=-M_landinggear_nose*9.80665
component_w[int((x_cg_mgear/l_fuselage)*n)]=-M_landinggear_main*9.80665
component_w[int((x_cg_wing_group/l_fuselage)*n)]=-M_wing_group*9.80665
   


lift_forces = [0] *n 
for i in range(int((X_wingbox_start/l_fuselage)*n),int((X_wingbox_end/l_fuselage)*n)):
    lift_forces[i]=Lift_mainwing/(int((X_wingbox_end/l_fuselage)*n)-int((X_wingbox_start/l_fuselage)*n))
    
lift_forces[int((x_cg_hwing/l_fuselage)*n)]=Lift_tail

forces_sum =  [0] *n 
for i in range(n):
    forces_sum[i]= payload_w[i]+fuselage_w[i]+component_w[i]+lift_forces[i]

V =  [0] *n 
for i in range(n):
    if i == 0:
        V[i]=0
    else:
        V[i]=V[i-1]+forces_sum[i-1]

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

plt.figure()
plt.plot(x,V)
plt.show()
