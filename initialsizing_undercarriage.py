# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:30:05 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt

# inputs
MTOW = 50000                                                                    # [kg] maximum take-off weight
x_cg = 10                                                                       # [m] x-location of the c.g
y_cg = 0                                                                        # [m] y-location of the c.g
z_cg = 3                                                                        # [m] z-location of the c.g
theta = 15                                                                      # [deg] scrape angle
beta = 15                                                                       # [deg] tip-back angle
phi = 5                                                                         # [deg] tip clearance angle
theta_rad = np.deg2rad(15)                                                      # [rad] scrape angle
beta_rad = np.deg2rad(15)                                                       # [rad] tip-back angle
phi_rad = np.deg2rad(5)                                                         # [rad] tip clearance angle
stroke = 0.2                                                                    # [m] shock absorber stroke



# constants
g = 9.80665                                                                     # [m/s^2] gravitational acceleration

N_mw = round(MTOW*g/120000)                                                     # [-] number of wheels mlg
N_nw = 2                                                                        # [-] number of wheels nlg
N_struts = 2                                                                    # [-] number of struts used

pressure_mw = 1.117E6                                                           # [Pa] tire pressure mlg
LCN = 46                                                                        # [-] Load classification number

weight_distribution = 0.08                                                      # [-] weight percentage on nose wheel

P_mw = (1-weight_distribution)*MTOW*g/N_mw                                      # [N] static loading on mw
P_nw = weight_distribution*MTOW*g/N_nw                                          # [N] static loading on nw

wheel_sizing_mw = (pressure_mw,P_mw/g)                                          # [-] data used to find dimentions mw

l_w = np.tan(beta_rad) * (z_cg - stroke)                                        # [m] mlg distance from c.g
l_n = (l_w*P_mw*N_mw)/(P_nw*N_nw)                                               # [m] nlg distance from c.g

x_w = x_cg + l_w                                                                # [m] x-location of mlg
x_n = x_cg - l_n                                                                # [m] x-location of nlg

plot = True
if plot:
    plt.plot([0,x_cg],[z_cg,z_cg])
    plt.axis('equal')
    plt.scatter(x_cg,z_cg, color = 'r')
    plt.scatter([x_n,x_w],[0,0])


