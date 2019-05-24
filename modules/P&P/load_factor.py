# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:16:35 2019

@author: Stijn, Rik
"""
import numpy as np
import matplotlib.pyplot as plt
from inputs.constants import *
from inputs.concept_1 import *

# inputs:
C_L_max = 1.5                                                                   # [-]
C_L_min = 1.5                                                                   # [-]
C_D_C_L_max = 0.05                                                              # [-]
MTOW = 65000                                                                    # [kg]
S = max(S)

def calc_V_S(MTOW,C_L_max,C_D_C_L_max,S,rho):
    W = MTOW *9.80665                                                           # [N]
    C_N_max = np.sqrt(C_L_max**2 + C_D_C_L_max**2)                              # [-]
    V_S = np.sqrt(2*W/(S * rho * C_N_max))                                      # [m/s]
    return V_S

def calc_n_lim_pos(MTOW):
    W = MTOW * 2.2046226218488                                                  # [lbs]
    n_lim_pos = 2.1 + 24000 / (W + 10000)                                       # [-]
    return n_lim_pos

def calc_V_C(W_S):
    W_S_psf = W_S / 47.88026                                                    # [psf]
    k_c = -0.0055*W_S_psf + 34.1                                                # [psf]
    return k_c * np.sqrt(W_S_psf)                                               # [m/s]

def calc_n_lim_neg(n_lim_pos):
    return 0.4* n_lim_pos                                                       # [-]

def calc_V_D(V_C):
    return 1.25 * V_C                                                           # [m/s]

def calc_V_A(V_S,V_C,n_lim):
    V_A = V_S * np.sqrt( n_lim)                                                 # [m/s]
    assert V_A < V_C
    return V_A

def calc_V_S_neg(MTOW,C_L_min,S,rho):
    W = MTOW *9.80665                                                           # [N]
    C_N_max_neg = 1.1 * C_L_min                                                 # [-]
    V_S_neg = np.sqrt(2*W/(S * rho * C_N_max_neg))                              # [m/s]
    return V_S_neg

V_S = calc_V_S(MTOW,C_L_max,C_D_C_L_max,S,rho)
n_lim_pos = calc_n_lim_pos(MTOW)
V_C = calc_V_C(W_S[0])
n_lim_neg = calc_n_lim_neg(n_lim_pos)
V_D = calc_V_D(V_C)
V_A = calc_V_A(V_S,V_C,n_lim_pos)
V_S_neg = calc_V_S_neg(MTOW,C_L_min,S,rho)

points = np.array([[0,0],
                  [0,1],
                  [V_S,1],
                  [V_A,n_lim_pos],
                  [V_D, n_lim_pos]])

fig = plt.figure(figsize = (12,5))
ax = fig.add_subplot(111)
ax.scatter(points[:,0],points[:,1])
ax.plot()

# positive Cn_max polynomial fit
pos_curve_fn = np.poly1d(np.polyfit([0, V_S, V_A], [0, 1.0, n_lim_pos], 2))
pos_curve_x = np.arange(0, V_A)
pos_curve_y = pos_curve_fn(pos_curve_x)
ax.plot(pos_curve_x, pos_curve_y)
