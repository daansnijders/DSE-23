# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:50:25 2019

@author: thong

Wing Analysis using Concept 1 values
Analyse canard
Analyse flap extension
Analyse wing tip extension
"""
import numpy as np
import matplotlib.pyplot as plt
"""
Parameters
Configuration 1
"""

b1 = 38.593274          # Span [m]
S1 = 135.403711         # Surface Area [m^2]
Cr1 = 5.359098          # Root Chord [m]
Ct1 = 1.657860          # Tip Chord [m]
MTOM1 = 58722.6         # lightest Maximum takeoff Mass [kg]

AR1 = b1**2/S1
TR1 = Ct1/Cr1

MTOW1 = MTOM1*9.81      # [N]
K = MTOW1/S1

"""
Parameters
Configuration 3
"""
MTOM2 = 68264.27        # Heaviest Maximum takeoff mass [kg]

MTOW2 = MTOM2*9.81      # [N]
S2 = MTOW2/K            # [m^2]

delta_S = S2-S1         # [m^2]

"""
Concept 1 (Canard)
"""
S_canard = delta_S
b_canard = (S_canard*AR1)**0.5
TR_canard = TR1
Cr_canard = 2*S_canard/((1+TR_canard)*b_canard)
Ct_canard = Cr_canard*TR_canard

"""
Concept 2
"""
S_concept2 = S2
b_concept2 = b1
AR_concept2 = b_concept2**2/S2
Ct_concept2 = Ct1
Cr_concept2 = 2*S_concept2/b_concept2 - Ct_concept2
Tr_concept2 = Cr_concept2/Ct_concept2

"""
Concept 3
"""
S_concept3 = S2
Cr_concept3 = Cr1
b_concept3 = (S_concept3*AR1)**0.5
Ct_concept3 = 2*S_concept3/b_concept3-Cr_concept3
#
def calculate_moments(b,S,MTOW, Cr, Ct):
    q = MTOW/S
    N = 10000
    step_size = b/(2*N)
    M = 0
    L_lst = []
    y_lst = []
    for i in range(N):
        y = b/2-i*step_size
        Ct_1 = Cr-(Cr-Ct)/(b/2)*y
        Cr_1 = Cr-(Cr-Ct)/(b/2)*(y-step_size)
        L = 0.5*(Cr_1+Ct_1)*step_size*q
        y_lst.append(y)
        L_lst.append(L)
    
    
    reac_force = -sum(L_lst)
    reac_moment = -sum(a*b for a,b in zip(L_lst,y_lst))
    
    delta_F = reac_force+MTOW/2
    L_lst = L_lst[::-1]
    y_lst = y_lst[::-1]
    moment_lst = []
    moment = 0
    for i in range(N):
        moment += L_lst[i]*y_lst[i]
        M = reac_moment + moment
        moment_lst.append(M)
        
    return [y_lst,moment_lst, delta_F]

concept1 = calculate_moments(b1,S1,MTOW1,Cr1,Ct1)
concept1b = calculate_moments(b_canard,S_canard,MTOW2-MTOW1,Cr_canard,Ct_canard)
concept2 = calculate_moments(b_concept2,S_concept2,MTOW2,Cr_concept2,Ct_concept2)
concept3 = calculate_moments(b_concept3,S_concept3,MTOW2,Cr_concept3,Ct_concept3)    

plt.figure()
plt.plot(concept1[0],concept1[1],label = "Concept 1")
plt.plot(concept1b[0],concept1b[1],label = "Concept 1 Canard")
plt.plot(concept2[0],concept2[1], label = "Concept 2")
plt.plot(concept3[0],concept3[1], label = "Concept 3")
plt.legend()    