# -*- coding: utf-8 -*-
"""
Created on Fri May 31 17:27:42 2019
@author: thong

Calculate geometrical property of the cross-section at y_location or chord
of the wing.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
"""
Parameters
------------------------------------------------------------------------------
"""

# Spar thickness [mm]
t_spar = 0.005

# Skin thickness [mm]
t_skin = 0.005

# Front spar location/chord [-]
spar_front = 0.25

# Rear spar location/chord [-]
spar_rear = 0.6

# Chord [m] example
chord = 5

# Stringer area [m^2]
stringer_area_upper = 0#314/10**6
stringer_area_lower = 0#314/10**6

"""
Read airfoil data
Obtain coordinate at upper and lower side of airfoil
------------------------------------------------------------------------------
"""
root_foil = open("66615.txt",'r')
tip_foil = open("65615.txt",'r')

root_foil_lines = root_foil.readlines()
tip_foil_lines = tip_foil.readlines()

root_foil_coor = []
tip_foil_coor = []
for i in root_foil_lines:
    x = float((i[:10].replace(',','.')))
    z = float((i[11:].replace(',','.')))
    root_foil_coor.append([x,z])
for j in tip_foil_lines:
    x = float((j[:10].replace(',','.')))
    z = float((j[11:].replace(',','.')))
    tip_foil_coor.append([x,z])

root_foil_coor = np.array(root_foil_coor)
tip_foil_coor = np.array(tip_foil_coor)

root_foil_coor[:,0] = root_foil_coor[:,0]-min(root_foil_coor[:,0])
tip_foil_coor[:,0] = tip_foil_coor[:,0]-min(tip_foil_coor[:,0])

root_foil_upper = root_foil_coor[:175]
root_foil_lower = root_foil_coor[175:]

tip_foil_upper = tip_foil_coor[:175]
tip_foil_lower = tip_foil_coor[175:]

"""
Read excel data
Obtain number of stringers and position of ribs
------------------------------------------------------------------------------
"""
component = pd.read_excel('wing_box_design.xlsx')

"""
Obtain cross-section properties
wing box area,full area,fuel area,moment of inertia,mass at cross section,limit stress 
------------------------------------------------------------------------------
"""

def get_cross_sec_prop(chord,y_loc,span):
    front_spar = 0.25*chord
    rear_spar = 0.6*chord
    N = 500
    step_size = (rear_spar-front_spar)/N
    
    x_lst = []
    z_lst = []
    

    y_frac = y_loc/(span/2)
    upper_stringer = int(np.interp(y_frac,component['y_frac'],component['stringer_up']))
    lower_stringer = int(np.interp(y_frac,component['y_frac'],component['stringer_down']))
    stringer_interval_upper = (rear_spar-front_spar)/(upper_stringer+1)
    stringer_interval_lower = (rear_spar-front_spar)/(lower_stringer+1)
    
    if y_frac>0.4:
        foil_upper = tip_foil_upper*chord
        foil_lower = tip_foil_lower*chord
    else:
        foil_upper = root_foil_upper*chord
        foil_lower = root_foil_lower*chord
    foil_upper=foil_upper[np.argsort(foil_upper[:, 0])]
    #------------------------------------------------------------------------------
    upper_frontspar_z = np.interp(front_spar,foil_upper[:,0],foil_upper[:,1])
    lower_frontspar_z = np.interp(front_spar,foil_lower[:,0],foil_lower[:,1])
    upper_rearspar_z = np.interp(rear_spar,foil_upper[:,0],foil_upper[:,1])
    lower_rearspar_z = np.interp(rear_spar,foil_lower[:,0],foil_lower[:,1])
    #-----------------------------Obtain total area and neutral_axis---------------
    tot_stringer_area = upper_stringer*stringer_area_upper+lower_stringer*stringer_area_lower
    front_spar_area = (upper_frontspar_z-lower_frontspar_z)*t_spar
    rear_spar_area = (upper_rearspar_z-lower_rearspar_z)*t_spar
    
    front_spar_area_z = front_spar_area * (upper_frontspar_z+lower_frontspar_z)/2
    rear_spar_area_z = rear_spar_area * (upper_rearspar_z+lower_rearspar_z)/2
    
    skin_box_area = 0
    skin_box_area_z = 0
    enclosed_area = 0
    for i in range(N):
        x1 = front_spar+i*step_size-0.5*step_size
        x2 = front_spar+i*step_size+0.5*step_size
        z1_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z2_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z1_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        z2_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        
        skin_box_area += abs(x2-x1)*t_skin + abs(x2-x1)*t_skin
        skin_box_area_z += abs(x2-x1)*t_skin*(z2_up+z1_up)/2+\
                            abs(x2-x1)*t_skin*(z2_low+z1_low)/2
        enclosed_area += abs(((z2_up+z1_up)/2 - (z2_low+z1_low)/2)*(x2-x1))    
        
        
    tot_stringer_area_z = 0
    for j in range(upper_stringer):
        xi = stringer_interval_upper*(j+1)+front_spar
        zi = np.interp(xi,foil_upper[:,0],foil_upper[0:,1])
        tot_stringer_area_z += stringer_area_upper*zi
        x_lst.append(xi)
        z_lst.append(zi)
        
        
    for k in range(lower_stringer):
        xi = stringer_interval_lower*(k+1)+front_spar
        zi = np.interp(xi,foil_lower[:,0],foil_lower[:,1])
        tot_stringer_area_z += stringer_area_lower*zi
        x_lst.append(xi)
        z_lst.append(zi)
        
    total_area = skin_box_area+tot_stringer_area+front_spar_area+rear_spar_area
    total_area_z = skin_box_area_z+tot_stringer_area_z+front_spar_area_z+rear_spar_area_z
    fuel_area = enclosed_area - total_area
    neutral_z = total_area_z/total_area
    max_z = max(foil_upper[:,1])-neutral_z
    min_z = min(foil_lower[:,1])-neutral_z
    
    #-------------------------Calculate area moment of inertia-----------------
    front_spar_moi = (1/12)*((upper_frontspar_z-lower_frontspar_z)**3)*t_spar + \
                        front_spar_area*((upper_frontspar_z+lower_frontspar_z)/2-
                                         neutral_z)**2

    rear_spar_moi = (1/12)*((upper_rearspar_z-lower_rearspar_z)**3)*t_spar + \
                        rear_spar_area*((upper_rearspar_z+lower_rearspar_z)/2-
                                         neutral_z)**2

    skin_moi = 0
    for i in range(N):
        x1 = front_spar+i*step_size-0.5*step_size
        x2 = front_spar+i*step_size+0.5*step_size
        z1_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z2_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z1_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        z2_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        skin_moi += abs(x2-x1)*t_skin*((z2_up+z1_up)/2-neutral_z)**2
        skin_moi += abs(x2-x1)*t_skin*((z2_low+z1_low)/2-neutral_z)**2
    
    upper_stringer_moi = 0
    for j in range(upper_stringer):
        xi = stringer_interval_upper*(j+1)+front_spar
        zi = np.interp(xi,foil_upper[:,0],foil_upper[0:,1])
        upper_stringer_moi += stringer_area_upper*(zi-neutral_z)**2
    
    lower_stringer_moi = 0
    for k in range(lower_stringer):
        xi = stringer_interval_lower*(k+1)+front_spar
        zi = np.interp(xi,foil_lower[:,0],foil_lower[:,1])
        lower_stringer_moi += stringer_area_upper*(zi-neutral_z)**2

    moi = front_spar_moi+rear_spar_moi+skin_moi+upper_stringer_moi+lower_stringer_moi
    return [enclosed_area,total_area,fuel_area,neutral_z,min_z,max_z,moi]

stuff = get_cross_sec_prop(5,1,1)