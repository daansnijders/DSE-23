# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:30:05 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt

# inputs
MTOW = 50000                                                                    # [kg] maximum take-off weight
x_cg = 17.5                                                                     # [m] x-location of the c.g
y_cg = 0                                                                        # [m] y-location of the c.g
z_cg = 3                                                                        # [m] z-location of the c.g
y_eng = 15                                                                      # [m] y-location of the engine
z_eng = 1.5                                                                     # [m] ground clearance of the engine
x_tail = 23                                                                     # [m] start of tail
l_f = 35                                                                        # [m] length of the fuselage
theta = 15                                                                      # [deg] scrape angle
beta = 15                                                                       # [deg] tip-back angle
phi = 5                                                                         # [deg] tip clearance angle
psi = 55                                                                        # [deg] overturn angle
theta_rad = np.deg2rad(theta)                                                   # [rad] scrape angle
beta_rad = np.deg2rad(beta)                                                     # [rad] tip-back angle
phi_rad = np.deg2rad(phi)                                                       # [rad] tip clearance angle
psi_rad = np.deg2rad(psi)                                                       # [rad] overturn angle
stroke = 0.2                                                                    # [m] shock absorber stroke



# constants
g = 9.80665                                                                     # [m/s^2] gravitational acceleration

N_mw = round(MTOW*g/120000)                                                     # [-] number of wheels mlg
N_nw = 2                                                                        # [-] number of wheels nlg
N_struts = 2                                                                    # [-] number of struts used

pressure_mw = 1.117E6                                                           # [Pa] tire pressure mlg
LCN = 46                                                                        # [-] load classification number

weight_distribution = 0.08                                                      # [-] weight percentage on nose wheel

P_mw = (1-weight_distribution)*MTOW*g/N_mw                                      # [N] static loading on mw
P_nw = weight_distribution*MTOW*g/N_nw                                          # [N] static loading on nw

wheel_sizing_mw = (pressure_mw,P_mw/g)                                          # [-] data used to find dimentions mw



# making eq of tangent lines
b1 = z_cg - np.tan(beta_rad - np.pi/2 ) * x_cg
b2 = - np.tan(theta_rad) * x_tail

x_mlg = (b2-b1)/(np.tan(beta_rad-np.pi/2) - np.tan(theta_rad))                  # [m] x-location of mlg
z_mlg = np.tan(theta_rad)*x_mlg + b2                                            # [m] z-location of mlg

l_w = x_mlg - x_cg                                                              # [m] mlg distance from c.g
l_n = (l_w*P_mw*N_mw)/(P_nw*N_nw)                                               # [m] nlg distance from c.g

x_nlg = x_cg - l_n                                                              # [m] x-location of nlg
y_nlg = 0                                                                       # [m] y-location of nlg
z_nlg = z_mlg                                                                   # [m] z-location of nlg


y_mlg = (l_n + l_w)/(np.sqrt((l_n**2 + np.tan(psi_rad)**2)/(z_cg + z_mlg) - 1))

ze = (y_mlg - y_eng) * - np.tan(phi_rad)

plot = True
if plot:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(xlabel = 'x', ylabel = 'z')
    ax.plot([0,l_f],[z_cg,z_cg])
    ax.axis('equal')

    ax.scatter(x_cg,z_cg, color = 'r')
    ax.scatter([x_mlg,x_nlg],[z_mlg,z_nlg])


