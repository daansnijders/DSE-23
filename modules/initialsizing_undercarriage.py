# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:30:05 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# inputs
                               
#y_eng = 0.326*b/2                                                               # [m] y-location of the engine
#d_eng = 2.006                                                                   # [m] diameter of the engine
#z_eng = -d_eng/2                                                                # [m] z-location of lowest part of the engine
#x_tail = 22.919                                                                 # [m] start of tail

                                                            

def get_P_mw(MTOW,N_mw,weight_distribution):
    return [(1-weight_distribution)*MTOW[i]*9.80665/N_mw for i in range(3)]                 

def get_P_nw(MTOW,N_nw,weight_distribution):
    return [weight_distribution*MTOW[i]*9.80665/N_nw for i in range(3)]

def get_wheel_sizing_mw(tire_pressure,P_mw):
    return [(tire_pressure/100,P_mw[i]/9.80665) for i in range(3)]

def get_wheel_sizing_nw(tire_pressure,P_nw):
    return [(tire_pressure/100,P_nw[i]/9.80665) for i in range(3)]



## making eq of tangent lines
#b1 = z_cg - np.tan(beta_rad - np.pi/2 ) * x_cg
#b2 = - np.tan(theta_rad) * x_tail
#
#x_mlg = (b2-b1)/(np.tan(beta_rad-np.pi/2) - np.tan(theta_rad))                  # [m] x-location of mlg
#z_mlg = np.tan(theta_rad)*x_mlg + b2 + stroke                                   # [m] z-location of mlg
#
#l_w = x_mlg - x_cg                                                              # [m] mlg distance from c.g
#l_n = (l_w*P_mw*N_mw)/(P_nw*N_nw)                                               # [m] nlg distance from c.g
#
#x_nlg = x_cg - l_n                                                              # [m] x-location of nlg
#y_nlg = 0                                                                       # [m] y-location of nlg
#z_nlg = z_mlg                                                                   # [m] z-location of nlg
#
#z_t = b * np.tan(dihedral_rad)                                                  # [m] z-location of wingtip
#
#y_mlg = (l_n + l_w)/(np.sqrt((l_n**2 * np.tan(psi_rad)**2)/(z_cg + z_mlg) - 1)) # [m] y-location1 of mlg
#y_mlg1 = y_eng - (z_eng - z_mlg)/np.tan(phi_rad)                                # [m] y-location2 of mlg
#y_mlg2 = b/2 - (z_t-z_mlg)/np.tan(phi_rad)                                      # [m] y-location3 of mlg
#y_mlg = max([y_mlg,y_mlg1,y_mlg2])                                              # [m] y-location of mlg
#
#print(x_mlg,y_mlg,z_mlg)
#print(x_nlg,y_nlg,z_nlg)
#
#def test_scrape_angle():
#    print((x_tail - x_mlg)*np.tan(theta_rad))
#    print(z_mlg)
#
#test_scrape_angle()
#
## Plotting stuff if one is interested.
#def set_aspect_equal_3d(ax):
#    """Fix equal aspect bug for 3D plots."""
#
#    xlim = ax.get_xlim3d()
#    ylim = ax.get_ylim3d()
#    zlim = ax.get_zlim3d()
#
#    from numpy import mean
#    xmean = mean(xlim)
#    ymean = mean(ylim)
#    zmean = mean(zlim)
#
#    plot_radius = max([abs(lim - mean_)
#                       for lims, mean_ in ((xlim, xmean),
#                                           (ylim, ymean),
#                                           (zlim, zmean))
#                       for lim in lims])
#
#    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
#    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
#    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])
#
#plot = False
#if plot:
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')    
#    ax.set(xlabel = 'x',ylabel = 'y', zlabel = 'z')
#    ax.plot([0,l_f],[0,0],[z_cg,z_cg])
#    ax.plot([x_wing,x_wing],[-b/2,b/2],[z_cg,z_cg])
#    ax.set_aspect('equal')
#    ax.scatter([x_cg],[y_cg],[z_cg],color = 'r')
#    ax.scatter([x_nlg,x_mlg,x_mlg], [y_nlg,-y_mlg,y_mlg], [z_nlg,z_mlg,z_mlg])
#    set_aspect_equal_3d(ax) 
