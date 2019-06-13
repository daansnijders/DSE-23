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
stringer_area = 0.0002
skin_yield = 372000000



x_cg_wing_group = 16.0
X_wingbox_start = 15.0
X_wingbox_end = 17.5
M_payload = 12500.0
M_fuselage = 9000.0
M_fittings = 9000.0
M_verticaltail = 1000.0
M_landinggear_nose = 200.0
M_landinggear_main = 1500.0
M_wing_group = 12000.0
M_horizontal_tail = 1500.0
M_fuel = 15000.0
M_total = M_payload+M_fuselage+M_fittings+M_verticaltail+M_landinggear_nose+M_landinggear_main+M_wing_group+M_horizontal_tail+M_fuel
#Lift_tail = -10000.0
#Lift_mainwing = M_total *g *nu - Lift_tail
l_fuselage = 35.0
x_cg_hwing = 34.0
x_cg_vwing = 34.0
x_cg_ngear = 5.0
x_cg_mgear = 16.0
Xfirst = 6.0
Xlast = 26.0

#l_fuselage = config1_class2.l_f
#x_cg_hwing = c1.x_cg_tail
#x_cg_vwing = c1.x_cg_tail
#x_cg_ngear = c1.x_nlg
#x_cg_mgear = c1.x_mlg
#Xfirst = c1.l_cockpit
#Xlast = Xfirst+c1.l_cabin
#M_horizontal_tail = config1_class2.M_horizontaltail


n= 1000
step_size = float(l_fuselage)/n

x = np.linspace(0,l_fuselage,n)



radius = [0] *n
for i in range(n):
    radius[i]=3.685*0.5
    
    
p_diff = isa(2438/3.281)[1] - isa(37000/3.281)[1]
skin_thickness = (p_diff* max(radius))/(skin_yield/safety_factor)
skin_thickness = 0.001


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

Lift_mainwing = (-m+x_cg_hwing*sum(weights_sum))/(np.average([X_wingbox_start, X_wingbox_end])-x_cg_hwing)
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
    if i == 0:
        V[i]=0
    else:
        V[i]=V[i-1]+forces_sum[i-1]

M =  [0] *n 
for i in range(n):
    for j in range(i):
        M[i]+=forces_sum[j]*(x[i]-step_size*j)

p_diff = isa(2438/3.281)[1] - isa(37000/3.281)[1]


stringer_no = [0] *n
MOI = [0] *n
stress_pressure_long =  [0] *n
stress_bending_max = [0] *n
stress_long_max = [1000000000000] *n
for i in range(n):
    while stress_long_max[i] > (stringer_yield/fitting_factor):       
        stringer_no[i]+=1
        MOI[i] = stringer_no[i]*(radius[i])**2 *stringer_area * 0.5 #+ np.pi*radius[i]**3*skin_thickness
    
        stress_pressure_long[i]=(p_diff*np.pi*(radius[i])**2)/(stringer_area*stringer_no[i] + 2*np.pi*radius[i]*skin_thickness)    
    
        stress_bending_max[i] = M[i]*radius[i]/MOI[i]    
       
        stress_long_max[i] = (abs(stress_bending_max[i])+stress_pressure_long[i])* safety_factor
     
frame_spacing = 0.5  
no_of_frames = int(l_fuselage/ frame_spacing)
n_in_frame = int(1.0/no_of_frames* n)
k = 0
stringer_no_new= [1]* n
for j in range(no_of_frames):
    for i in range(k, k + n_in_frame):
        stringer_no_new[i] = max(stringer_no[k:k+n_in_frame])
    k+= n_in_frame
    
    
    
    
    
    

plt.figure()
plt.plot(x,stringer_no)
plt.show()

plt.figure()
plt.plot(x,stringer_no_new)
plt.show()


plt.figure()
plt.plot(x, V)
plt.show()

plt.figure()
plt.plot(x, M)
plt.show()

######################################BUCKILING LOAD FACTOR TEST#####################################    
nu = 1       
buckling_crit = 100000000000000.0 
start = True
while start:
    nu+=0.01
    
        
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
    
    Lift_mainwing = (-m+x_cg_hwing*sum(weights_sum))/(np.average([X_wingbox_start, X_wingbox_end])-x_cg_hwing)
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
        if i == 0:
            V[i]=0
        else:
            V[i]=V[i-1]+forces_sum[i-1]
    
    M =  [0] *n 
    for i in range(n):
        for j in range(i):
            M[i]+=forces_sum[j]*(x[i]-step_size*j)
    
    
    
    buckling_crit = [0]*n
    for i in range(n):      
        MOI[i] = stringer_no_new[i]*(radius[i])**2 *stringer_area * 0.5 + np.pi*radius[i]**3*skin_thickness
       
        Y_mod = 72 * 10**9
        v = 0.33
        b_skin =  2* np.pi /stringer_no_new[i] * radius[i]
        buckling_crit[i] = (skin_thickness/b_skin)**2 * (np.pi**2 * 7 * Y_mod)/(12*(1-v**2))
        
    
        stress_pressure_long[i]=(p_diff*np.pi*(radius[i])**2)/(stringer_area*stringer_no[i] + 2*np.pi*radius[i]*skin_thickness)    
    
        stress_bending_max[i] = M[i]*radius[i]/MOI[i]    
       
        stress_long_max[i] = (abs(stress_bending_max[i])-stress_pressure_long[i])* safety_factor
        
    for i in range(n):
        if stress_long_max[i] > buckling_crit[i]:
            start = False
        
print "Load factor:", nu        

    
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

##########Titanium Ti-6Al-4V aged##########
delta = 0.0001
b = 0.5 *t +0.25 *t + delta
Mom = sigma*A*(b/2)
bolt_d = D - delta*2
bolt_rupture = 1020000000
I_bolt = np.pi/4 * (bolt_d/2)**4
MS = ((Mom*bolt_d)/(2*I_bolt))/bolt_rupture  -1
print "MS:", MS


#################################Frame spacing################################
width = 0.0625
r_frame = max(radius)
I_frame = np.pi/4 * (r_frame**4 - (r_frame-width)**4)
######frame: Aluminum 2024-T81########
Y_mod = 72 * 10**9
v = 0.33
frame_yield = 372000000
L = abs(min(M))*(r_frame*2)**2/(16000*Y_mod*I_frame)
b_skin =  2* np.pi /180 *r_frame
k_c = (b_skin/0.001)**2 * (max(stress_long_max) * 12 * (1-v**2))/(Y_mod * np.pi**2)
print (k_c)





