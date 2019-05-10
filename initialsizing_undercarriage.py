# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:30:05 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# inputs
MTOW = 48229.58763                                                              # [kg] maximum take-off weight
x_cg = 15.95274574                                                              # [m] x-location of the c.g
y_cg = 0                                                                        # [m] y-location of the c.g
z_cg = 3.685/2                                                                  # [m] z-location of the c.g
b = 34.9756475                                                                  # [m] wingspan                                                            
y_eng = 0.326*b/2                                                               # [m] y-location of the engine
d_eng = 2.006                                                                   # [m] diameter of the engine
z_eng = -d_eng/2                                                                # [m] z-location of lowest part of the engine
x_tail = 22.919                                                                 # [m] start of tail
l_f = 35.815                                                                    # [m] length of the fuselage
theta = 15                                                                      # [deg] scrape angle
beta = 17                                                                       # [deg] tip-back angle
phi = 5                                                                         # [deg] tip clearance angle
psi = 55                                                                        # [deg] overturn angle
theta_rad = np.deg2rad(theta)                                                   # [rad] scrape angle
beta_rad = np.deg2rad(beta)                                                     # [rad] tip-back angle
phi_rad = np.deg2rad(phi)                                                       # [rad] tip clearance angle
psi_rad = np.deg2rad(psi)                                                       # [rad] overturn angle
stroke = 0.3                                                                    # [m] shock absorber stroke
dihedral = 2.403196373                                                          # [deg] dihedral angle
dihedral_rad = np.deg2rad(dihedral)                                             # [rad] dihedral angle
x_wing = 15


# constants
g = 9.80665                                                                     # [m/s^2] gravitational acceleration

N_mw = 4                                                                        # [-] number of wheels mlg
N_nw = 2                                                                        # [-] number of wheels nlg
N_struts = 2                                                                    # [-] number of struts used

LCN = 45                                                                        # [-] load classification number
tire_pressure = 430 * np.log(LCN) - 680                                         # [Pa] tire pressure mlg

weight_distribution = 0.08                                                      # [-] weight percentage on nose wheel

P_mw = (1-weight_distribution)*MTOW*g/N_mw                                      # [N] static loading on mw
P_nw = weight_distribution*MTOW*g/N_nw                                          # [N] static loading on nw


wheel_sizing_mw = (tire_pressure/100,P_mw/g)                                    # [-] data used to find dimentions mw
wheel_sizing_nw = (tire_pressure/100,P_nw/g)                                    # [-] data used to find dimentions nw



# making eq of tangent lines
b1 = z_cg - np.tan(beta_rad - np.pi/2 ) * x_cg
b2 = - np.tan(theta_rad) * x_tail

x_mlg = (b2-b1)/(np.tan(beta_rad-np.pi/2) - np.tan(theta_rad))                  # [m] x-location of mlg
z_mlg = np.tan(theta_rad)*x_mlg + b2 + stroke                                   # [m] z-location of mlg

l_w = x_mlg - x_cg                                                              # [m] mlg distance from c.g
l_n = (l_w*P_mw*N_mw)/(P_nw*N_nw)                                               # [m] nlg distance from c.g

x_nlg = x_cg - l_n                                                              # [m] x-location of nlg
y_nlg = 0                                                                       # [m] y-location of nlg
z_nlg = z_mlg                                                                   # [m] z-location of nlg

z_t = b * np.tan(dihedral_rad)                                                  # [m] z-location of wingtip

y_mlg = (l_n + l_w)/(np.sqrt((l_n**2 * np.tan(psi_rad)**2)/(z_cg + z_mlg) - 1)) # [m] y-location1 of mlg
y_mlg1 = y_eng - (z_eng - z_mlg)/np.tan(phi_rad)                                # [m] y-location2 of mlg
y_mlg2 = b/2 - (z_t-z_mlg)/np.tan(phi_rad)                                      # [m] y-location3 of mlg
y_mlg = max([y_mlg,y_mlg1,y_mlg2])                                              # [m] y-location of mlg

print(x_mlg,y_mlg,z_mlg)
print(x_nlg,y_nlg,z_nlg)

def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])

plot = True
if plot:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')    
    ax.set(xlabel = 'x',ylabel = 'y', zlabel = 'z')
    ax.plot([0,l_f],[0,0],[z_cg,z_cg])
    ax.plot([x_wing,x_wing],[-b/2,b/2],[z_cg,z_cg])
    ax.set_aspect('equal')
    ax.scatter([x_cg],[y_cg],[z_cg],color = 'r')
    ax.scatter([x_nlg,x_mlg,x_mlg], [y_nlg,-y_mlg,y_mlg], [z_nlg,z_mlg,z_mlg])
    set_aspect_equal_3d(ax) 
