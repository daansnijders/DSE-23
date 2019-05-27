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
S = max(S)                                                                      # [m^2]

def calc_n_lim_pos(MTOW):
    MTOW_lbf = MTOW * 2.2046226218488                                           # [lbf]
    n_lim_pos = 2.1 + 24000 / (MTOW_lbf + 10000)                                # [-]
    return n_lim_pos

def calc_n_lim_neg(n_lim_pos):
    n_lim_neg = max([0.4 * n_lim_pos, 1])                                       # [-]
    return n_lim_neg                                                            # [-]

def calc_V_S(W_S,C_L_max,C_D_C_L_max,rho):
    C_N_max = np.sqrt(C_L_max**2 + C_D_C_L_max**2)                              # [-]
    V_S = np.sqrt(2* W_S / (rho * C_N_max))                                     # [m/s]
    return V_S

def calc_V_C(W_S):
    W_S_psf = W_S / 47.88026                                                    # [psf]
    k_c = -0.0055*W_S_psf + 34.1                                                
    V_C = k_c * np.sqrt(W_S_psf) * kts_to_ms                                    # [m/s]
    return max([V_cruise,V_C])                                                  # [m/s]

def calc_V_D(V_C):
    return 1.25 * V_C                                                           # [m/s]

def calc_V_A(V_S,V_C,n_lim):
    V_A = V_S * np.sqrt(n_lim)                                                  # [m/s]
    assert V_A < V_C
    return V_A

def calc_V_S_neg(W_S,C_L_min,rho):
    W_S_psf = W_S / 47.88026                                                    # [psf]
    rho_scf = rho * 0.00194032                                                  # [sl/ft^3]
    C_N_max_neg = 1.1 * C_L_min                                                 # [-]
    V_S_neg = np.sqrt(2*W_S_psf/(rho_scf * C_N_max_neg)) * kts_to_ms            # [m/s]
    return V_S_neg

def calc_gust_C(W_S,MAC,rho,H_ft,V_C,C_L_alpha,g):
    W_S_psf = W_S / 47.88026                                                    # [psf]
    rho_sl = rho * 0.0019403203                                                 # [slug/ft^3]
    MAC_ft = MAC * m_to_ft                                                      # [ft]
    g_ft = g * 3.28084                                                          # [ft/s^2]
    V_keas = V_C / kts_to_ms                                                    # [KEAS]

    mu_g = 2 * W_S_psf / (rho_sl * MAC_ft * g_ft * C_L_alpha)                   # [-]
    K_g = 0.88 * mu_g / (5.3 + mu_g)                                            # [-]
    U_de_C = 66.67 - 0.000833 * H_ft                                            # [ft/s]
    n_lim_C = 1 + (K_g * U_de_C * V_keas * C_L_alpha)/(498 * W_S_psf)           # [-]
    return n_lim_C

def calc_gust_D(W_S,MAC,rho,H_ft,V_D,C_L_alpha,g):
    W_S_psf = W_S / 47.88026                                                    # [psf]
    rho_sl = rho * 0.0019403203                                                 # [slug/ft^3]
    MAC_ft = MAC * m_to_ft                                                      # [ft]
    g_ft = g * 3.28084                                                          # [ft/s^2]
    V_keas = V_D / kts_to_ms                                                    # [KEAS]

    mu_g = 2 * W_S_psf / (rho_sl * MAC_ft * g_ft * C_L_alpha)                   # [-]
    K_g = 0.88 * mu_g / (5.3 + mu_g)                                            # [-]
    U_de_D = 33.34 - 0.000417 * H_ft                                            # [ft/s]
    n_lim_D = 1 + (K_g * U_de_D * V_keas * C_L_alpha)/(498 * W_S_psf)           # [-]
    return n_lim_D
    

V_S = calc_V_S(W_S[0],C_L_max,C_D_C_L_max,rho)                                  # [m/s]
n_lim_pos = calc_n_lim_pos(MTOW)                                                # [-]
V_C = calc_V_C(W_S[0])                                                          # [m/s]
n_lim_neg = calc_n_lim_neg(n_lim_pos)                                           # [-]
V_D = calc_V_D(V_C)                                                             # [m/s]
V_A = calc_V_A(V_S,V_C,n_lim_pos)                                               # [m/s]
V_S_neg = calc_V_S_neg(MTOW,C_L_min,rho)                                        # [m/s]
n_lim_C = calc_gust_C(W_S[0],MAC[2],rho,H_ft,V_C,CL_alpha[0],g)                 # [-]
n_lim_D = calc_gust_D(W_S[0],MAC[2],rho,H_ft,V_D,CL_alpha[0],g)                 # [-]

# tests
assert V_S < V_A < V_C < V_D

# positive Cn_max polynomial fit
pos_curve_fn = np.poly1d(np.polyfit([0, V_S, V_A], [0, 1.0, n_lim_pos], 2))
pos_curve_x1 = np.linspace(0, V_S,50)
pos_curve_y1 = pos_curve_fn(pos_curve_x1)
pos_curve_x2 = np.linspace(V_S, V_A,50)
pos_curve_y2 = pos_curve_fn(pos_curve_x2)

# negative Cn_max polynomial fit
if n_lim_neg > 1:
    neg_curve_x1 = np.linspace(0, V_S_neg, 50)
    neg_curve_x2 = np.linspace(V_S_neg, V_A, 50)
    neg_curve_fn = np.poly1d(np.polyfit([0, V_S_neg, V_A], [0, -1.0, -n_lim_neg], 2))
    neg_curve_y1 = neg_curve_fn(neg_curve_x1)
    neg_curve_y2 = neg_curve_fn(neg_curve_x2)
else:
    neg_curve_x1 = np.linspace(0, V_S, 50)
    neg_curve_x2 = np.linspace(V_S, V_A, 50)
    neg_curve_y1 = pos_curve_fn(neg_curve_x1) *-1
    neg_curve_y2 = np.ones(len(neg_curve_x2)) * -1

# line 0,1 to Vc
vc_gust_line_fn = np.polyfit([0, V_C], [1, n_lim_C], 1)
vc_gust_line_intersect = (n_lim_pos - 1)/vc_gust_line_fn[0]
vc_gust_line_x_1 = np.linspace(0, vc_gust_line_intersect, 50)
vc_gust_line_y_1 = np.linspace(1, n_lim_pos, 50)
vc_gust_line_x_2 = np.linspace(vc_gust_line_intersect, V_C, 50)
vc_gust_line_y_2 = np.linspace(n_lim_pos, n_lim_C, 50)

# line Vc,NlimC to Vd,NlimD
vc_vd = np.polyfit([V_C, V_D], [n_lim_C, n_lim_D], 1)
# print(vc_vd)

points = np.array([[0,0],
                  [0,1],
                  [V_S,1],
                  [V_A,n_lim_pos],
                  [V_D, n_lim_pos],
                  [V_S,-n_lim_neg],
                   [V_C, n_lim_C],
                   [V_D, n_lim_D]])

fig = plt.figure(figsize = (12,5))
ax = fig.add_subplot(111)
ax.scatter(points[:,0],points[:,1])

ax.plot(neg_curve_x1, neg_curve_y1, color = 'C0', linestyle = ':')
ax.plot(neg_curve_x2, neg_curve_y2, color = 'C0')

ax.plot(pos_curve_x1, pos_curve_y1, color = 'C0', linestyle = ':')
ax.plot(pos_curve_x2, pos_curve_y2, color = 'C0')

ax.plot(vc_gust_line_x_1, vc_gust_line_y_1, color = 'C0', linestyle = ':')
ax.plot(vc_gust_line_x_2, vc_gust_line_y_2, color = 'C0')
