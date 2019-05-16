# -*- coding: utf-8 -*-
"""
Created on Wed May 15 13:44:47 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt

b = 42.32651161
S = 162.8666896

Cr = 5.877498919
Ct = 1.818230465
L_total = 35000*9.80665*2.1
lambda_2_deg = 23.70581453 

sigma_ult = 159E6                                                               # [Pa]
thickness = 0.2                                                                 # [m]
M_eng = 2177                                                                    # [kg]
y_eng = 6.3                                                                     # [m]

def calc_chord_length(y, b=b, Cr=Cr, Ct=Ct):
    a = (Ct-Cr) / (b/2)
    return a*y+Cr

def calc_area(length,C_l, C_s):
    A = C_s * length + 0.5*(C_l - C_s)*length
    return A

dy = 200
lift_distr = np.linspace(0,b/2,dy,retstep=True)                                 # [m]
A = []
for length in lift_distr[0][:-1]:
    C_s = calc_chord_length(length + lift_distr[1])
    C_l = calc_chord_length(length)
    area = calc_area(lift_distr[1],C_l,C_s)
    A.append([area, length + 0.5*lift_distr[1]])
    
F_total = []
for i in range(len(A)):
    if A[i][1] > y_eng and not [y_eng,M_eng * -9.80665] in F_total:
        F_total.append([y_eng,M_eng * -9.80665])
    F_total.append([A[i][1],A[i][0]/(S/2)*L_total])
    
F_total = np.array(F_total)[::-1]


# Moment calculation
def calc_shear(y,F_total):
    V = 0
    i=0
    while not (i+1) > len(F_total) and y < F_total[i][0]:
        i +=1
    for forces in range(i):
        V += F_total[forces][1]   
    return V

def calc_moment(y,F_total):
    M = 0
    i = 0
    while not (i+1) > len(F_total) and y < F_total[i][0]:
        i +=1
    for forces in range(i):
        M += (F_total[forces][0] - y) * F_total[forces][1]    
    return [M,y]

bending_distr = np.linspace(0,b/2,200)

M = []
for y in bending_distr:
    M.append(calc_moment(y,F_total))
M = np.array(M)





