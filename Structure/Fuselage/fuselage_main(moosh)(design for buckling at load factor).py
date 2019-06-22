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
#import sys
#sys.path.insert(0,'/Users/thong/Documents/TU Delft DSE/DSE-23')
#import inputs.concept_1 as c1
#import modules.CG.CG_func as cgfunc
#import modules.CG.class2_CG as cl2cg
#
#from modules.class2 import *
#from main import *


fitting_factor = 1.15
safety_factor = 1.5
nu = 3
g = 9.80665
######stringers: Aluminum 2024-T81########
stringer_yield = 372000000
stringer_comp = 262000000
stringer_area = 0.0002
skin_yield = 627*10**6



X_wingbox_start = 11.12#
X_wingbox_end = 13.84#
x_cg_wing_group = np.average([X_wingbox_start, X_wingbox_end])#
M_payload = 10306.40#
M_fuselage = 7665.99#
M_fittings = 7934.3#
M_verticaltail = 58.97#
M_landinggear_nose = 364.07#
M_landinggear_main = 1967.13#
M_wing_group = 17076.11#
M_horizontal_tail = 75.73#
M_fuel = 6742.59#
M_total = M_payload+M_fuselage+M_fittings+M_verticaltail+M_landinggear_nose+M_landinggear_main+M_wing_group+M_horizontal_tail+M_fuel
l_fuselage = 30.13#
x_cg_hwing = 29.51#
x_cg_vwing = 29.84#
x_cg_ngear = 2.0#
x_cg_mgear = 15.34#
Xfirst = 7.09#
Xlast = 20.91#
wing_moment = 214599.7501376995/(3.0)




n= 2000
step_size = float(l_fuselage)/n

x = np.linspace(0,l_fuselage,n)



radius = [0] *n
for i in range(n):
    radius[i]=3.685*0.5
    
    
p_diff = isa(2438/3.281)[1] - isa(37000/3.281)[1]
skin_thickness = (p_diff* max(radius))/(skin_yield/safety_factor)
print skin_thickness
skin_thickness = 0.0015


payload_w = [0] *n
for i in range(int((Xfirst/l_fuselage)*n),int((Xlast/l_fuselage)*n)):
    payload_w[i]=-M_payload* g*nu/(int((Xlast/l_fuselage)*n)-int((Xfirst/l_fuselage)*n)) 
    
fuselage_w = [0] *n
for i in range(n):
    fuselage_w[i]=-(M_fuselage+M_fittings)*g*nu/n
    
component_w = [0] *n 
component_w[int((x_cg_hwing/l_fuselage)*n)]=-M_horizontal_tail* g*nu
component_w[int((x_cg_vwing/l_fuselage)*n)]=-M_verticaltail* g*nu
component_w[int((x_cg_ngear/l_fuselage)*n)]=-M_landinggear_nose* g*nu
component_w[int((x_cg_mgear/l_fuselage)*n)]=-M_landinggear_main* g*nu
#component_w[int((x_cg_wing_group/l_fuselage)*n)]=-(M_wing_group+ 1.00* M_fuel)* g*nu
for i in range(int((X_wingbox_start/l_fuselage)*n),int((X_wingbox_end/l_fuselage)*n)):
    component_w[i]=-(M_wing_group+ 1.00* M_fuel)* g*nu/(int((X_wingbox_end/l_fuselage)*n)-int((X_wingbox_start/l_fuselage)*n))   

weights_sum =  [0] *n 
for i in range(n):
    weights_sum[i]= payload_w[i]+fuselage_w[i]+component_w[i]

#m =  [0] *n 
#for i in range(n):
#    for j in range(i):
#        m[i]+=weights_sum[j]*(x[i]-step_size*j)
#        
m = 0     
for i in range(n):
    m+=weights_sum[i]*x[i]

Lift_mainwing = (-(m-wing_moment*nu)+x_cg_hwing*sum(weights_sum))/(np.average([X_wingbox_start, X_wingbox_end])-x_cg_hwing)
Lift_tail = -sum(weights_sum)-Lift_mainwing

lift_forces = [0] *n 
for i in range(int((X_wingbox_start/l_fuselage)*n),int((X_wingbox_end/l_fuselage)*n)):
    lift_forces[i]=Lift_mainwing/(int((X_wingbox_end/l_fuselage)*n)-int((X_wingbox_start/l_fuselage)*n))
    
lift_forces[int((x_cg_hwing/l_fuselage)*n)]=Lift_tail

forces_sum =  [0] *n 
for i in range(n):
    forces_sum[i]= payload_w[i]+fuselage_w[i]+component_w[i]+lift_forces[i]

V =  [0] *n 
for i in range(n):
    V[i]=V[i-1]+forces_sum[i]


M =  [0] *n 
for i in range(n):
    if i > int((x_cg_wing_group/l_fuselage)*n):
        for j in range(i):
            M[i]+=forces_sum[j]*(x[i]-step_size*j)
        M[i] += wing_moment*nu        
    else:
        for j in range(i):
            M[i]+=forces_sum[j]*(x[i]-step_size*j)
#M[int((x_cg_wing_group/l_fuselage)*n)]+=wing_moment*nu

p_diff = isa(2438/3.281)[1] - isa(37000/3.281)[1]


stringer_no = [0] *n
MOI = [0] *n
stress_pressure_long =  [0] *n
stress_bending_max = [0] *n
stress_long_max = [1000000000000] *n
critical = 0
for i in range(n):
    while stress_long_max[i] > (stringer_yield):       
        stringer_no[i]+=1
        MOI[i] = stringer_no[i]*(radius[i])**2 *stringer_area * 0.5 #+ np.pi*radius[i]**3*skin_thickness
       
    
        stress_pressure_long[i]=(p_diff*np.pi*(radius[i])**2)/(stringer_area*stringer_no[i] + 2*np.pi*radius[i]*skin_thickness)    
    
        stress_bending_max[i] = M[i]*radius[i]/MOI[i]    
       
        stress_long_max[i] = (abs(stress_bending_max[i])+stress_pressure_long[i])* safety_factor





plt.figure()
plt.plot(x, V)
plt.xlabel("x location in fuselage [m]")
plt.ylabel("Internal shear [N]")
plt.title("Internal shear across fuselage length")
plt.show()

plt.figure()
plt.plot(x, M)
plt.xlabel("x location in fuselage [m]")
plt.ylabel("Internal bending moment[Nm]")
plt.title("Internal bending moment across fuselage length")
plt.show()
     
#######################################BUCKLING CALCULATION: REQUIRED NO OF STRINGERS###################################
nu = 2
payload_w = [0] *n
for i in range(int((Xfirst/l_fuselage)*n),int((Xlast/l_fuselage)*n)):
    payload_w[i]=-M_payload* g*nu/(int((Xlast/l_fuselage)*n)-int((Xfirst/l_fuselage)*n)) 
    
fuselage_w = [0] *n
for i in range(n):
    fuselage_w[i]=-(M_fuselage+M_fittings)*g*nu/n
    
component_w = [0] *n 
component_w[int((x_cg_hwing/l_fuselage)*n)]=-M_horizontal_tail* g*nu
component_w[int((x_cg_vwing/l_fuselage)*n)]=-M_verticaltail* g*nu
component_w[int((x_cg_ngear/l_fuselage)*n)]=-M_landinggear_nose* g*nu
component_w[int((x_cg_mgear/l_fuselage)*n)]=-M_landinggear_main* g*nu
#component_w[int((x_cg_wing_group/l_fuselage)*n)]=-(M_wing_group+ 1.00* M_fuel)* g*nu
for i in range(int((X_wingbox_start/l_fuselage)*n),int((X_wingbox_end/l_fuselage)*n)):
    component_w[i]=-(M_wing_group+ 1.00* M_fuel)* g*nu/(int((X_wingbox_end/l_fuselage)*n)-int((X_wingbox_start/l_fuselage)*n))   

weights_sum =  [0] *n 
for i in range(n):
    weights_sum[i]= payload_w[i]+fuselage_w[i]+component_w[i]

#m =  [0] *n 
#for i in range(n):
#    for j in range(i):
#        m[i]+=weights_sum[j]*(x[i]-step_size*j)
#        
m = 0     
for i in range(n):
    m+=weights_sum[i]*x[i]

Lift_mainwing = (-(m-wing_moment*nu)+x_cg_hwing*sum(weights_sum))/(np.average([X_wingbox_start, X_wingbox_end])-x_cg_hwing)
Lift_tail = -sum(weights_sum)-Lift_mainwing

lift_forces = [0] *n 
for i in range(int((X_wingbox_start/l_fuselage)*n),int((X_wingbox_end/l_fuselage)*n)):
    lift_forces[i]=Lift_mainwing/(int((X_wingbox_end/l_fuselage)*n)-int((X_wingbox_start/l_fuselage)*n))
    
lift_forces[int((x_cg_hwing/l_fuselage)*n)]=Lift_tail

forces_sum =  [0] *n 
for i in range(n):
    forces_sum[i]= payload_w[i]+fuselage_w[i]+component_w[i]+lift_forces[i]

V =  [0] *n 
for i in range(n):
    V[i]=V[i-1]+forces_sum[i]

M =  [0] *n 
for i in range(n):
    if i > int((x_cg_wing_group/l_fuselage)*n):
        for j in range(i):
            M[i]+=forces_sum[j]*(x[i]-step_size*j)
        M[i] += wing_moment*nu        
    else:
        for j in range(i):
            M[i]+=forces_sum[j]*(x[i]-step_size*j)


p_diff = isa(2438/3.281)[1] - isa(37000/3.281)[1]


stringer_no_buckle = [0] *n
MOI = [0] *n
stress_pressure_long =  [0] *n
stress_bending_max = [0] *n
stress_long_max = [1000000000000] *n
critical = 0
for i in range(n):
    while stress_long_max[i] > critical:       
        stringer_no_buckle[i]+=1
        MOI[i] = stringer_no_buckle[i]*(radius[i])**2 *stringer_area * 0.5 + np.pi*radius[i]**3*skin_thickness
       
        Y_mod = 62.5 * 10**9
        v = 0.058
        v = 0.33
        b_skin =  2* np.pi /stringer_no_buckle[i] * radius[i]
        buckling_crit = (skin_thickness/b_skin)**2 * (np.pi**2 * 7 * Y_mod)/(12*(1-v**2))
        critical = min([buckling_crit, stringer_comp])
    
        stress_pressure_long[i]=(p_diff*np.pi*(radius[i])**2)/(stringer_area*stringer_no_buckle[i] + 2*np.pi*radius[i]*skin_thickness)    
    
        stress_bending_max[i] = M[i]*radius[i]/MOI[i]    
       
        stress_long_max[i] = (abs(stress_bending_max[i])-stress_pressure_long[i])* safety_factor
     











frame_spacing = 0.5  
no_of_frames = int(l_fuselage/ frame_spacing)
n_in_frame = int(1.0/no_of_frames* n)
k = 0
stringer_no_new1= [1]* n
for j in range(no_of_frames):
    for i in range(k, k + n_in_frame):
        stringer_no_new1[i] = max([max(stringer_no[k:k+n_in_frame]),max(stringer_no_buckle[k:k+n_in_frame])])
    k+= n_in_frame
        
    
######################################JOINT CALCULATION########################
radius_joint = 3.685*0.5
x_cut = 6.0
joint_no = 12.0
t = 0.0001  
sigma = 1000000000000000.0
sigma_critical = 0.0
while sigma > sigma_critical:
    t+=0.0001
    D = 4 * t  
    a = 1.5 * D
    w = 2 * a
    A = w * t
    ###########Titanium Ti-6Al-4V aged##########
    sigma_ult = 1100000000.0    
    sigma_y = 1020000000.0
    sigma_yt = 1100000000.0
    Kbr = 1.1
    Kt = 0.94
    sigma_cr1 = Kbr * sigma_ult * D/w
    sigma_cr2 = Kt * (w-D) * sigma_ult/w
    factor = min(sigma_cr1, sigma_cr2) *w*t/(D*t* sigma_ult)
    if factor <= 3.0 and factor >= 1.0:
        C = -0.2105* factor + 1.3473
    else:
        print ("ERROR with C value")
    sigma_cr3 = C * (sigma_y/sigma_ult) * min(sigma_cr1, sigma_cr2)
    sigma_critical = (min(sigma_cr1, sigma_cr2, sigma_cr3))/fitting_factor
    
    MOI_joint = joint_no* radius_joint**2 *0.5 * A
    sigma_press=(p_diff*np.pi*(radius_joint)**2)/(A*joint_no)  
    sigma_bend = (M[int((x_cut/l_fuselage)*n)]* radius_joint/MOI_joint)* safety_factor
    sigma = sigma_press + abs(sigma_bend)
    
    

print t


skin_density = 1540
stringer_density = 2780
skin_weight = l_fuselage * max(radius) *2 *np.pi * skin_density * skin_thickness
stringer_weight = np.trapz(stringer_no_new1, x)* stringer_area * stringer_density
    


plt.figure()
plt.plot(x,stringer_no_new1)
plt.show()

plt.figure()
plt.plot(x,stringer_no)
plt.show()

plt.figure()
plt.plot(x,stringer_no_buckle)
plt.show()


