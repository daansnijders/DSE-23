# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:16:32 2019

@author: moosh



Assumptions:
Three modes of failure:
    a)net section tension;
    b)shear tearout, assuming all the load to be transmitted on "40 degree planes"
    c)bearing
"""

import numpy as np
#import stress from fuselage_main
#
radius = 3.685
joint_no = 12
t = 0.01     
D = 4 * t  
a = 2 * D
w = 2 * a
A = w * t
sigma_ult = 1100000000
sigma_y = 1020000000
sigma_yt = 1100000000
Kbr = 1.1
Kt = 0.9
sigma_cr1 = Kbr * sigma_ult * D/w
sigma_cr2 = Kt * (w-D) * sigma_ult/w
factor = min(sigma_cr1, sigma_cr2) *w*t/(D*t* sigma_ult)
if factor <= 3.0 and factor >= 1.0:
    C = -0.2105* factor + 1.3473
else:
    print ("ERROR with C value")
sigma_cr3 = C * (sigma_y/sigma_ult) * min(sigma_cr1, sigma_cr2)
sigma_critical = min(sigma_cr1, sigma_cr2, sigma_cr3)

MOI_joint = joint_no* radius**2 *0.5 * A
sigma = M[i]* radius[i]/MOI_joint
if sigma > sigma_crit:
    print("we're fucked")

