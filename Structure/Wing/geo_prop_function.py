# -*- coding: utf-8 -*-
"""
Created on Fri May 31 17:27:42 2019
@author: thong

Calculate geometrical property of the cross-section at y_location or chord
of the wing.
"""
import numpy as np
import pandas as pd
"""
Parameters
------------------------------------------------------------------------------
"""

# Spar thickness [m]
t_spar = 0.008

# Skin thickness [m]
t_skin = 0.008

# Front spar location/chord [-]
spar_front = 0.15#0.1

# Rear spar location/chord [-]
spar_rear = 0.6

# Chord [m] example
chord = 5

stringer_area_upper = 600/10**6
stringer_area_lower = 600/10**6


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
component = pd.read_excel('wing_box_design.xlsx','design')

ks = pd.read_excel('wing_box_design.xlsx','ks_hinged')
ribs = component['ribs'] > 0
ribs_loc = component[ribs]['y_frac']
ribs_loc = np.array(ribs_loc)

# Materials properties for web [Pa]
E_web = 71*10**9
shear_strength_web = 310*10**6
web_poisson_ratio = 0.33

"""
Obtain web buckling stress
------------------------------------------------------------------------------
"""

def get_web_buckling(y_loc,spar_height,span):

    y_frac = y_loc/(span/2)
    i = 0
    for i in range(1,len(ribs_loc)):
        if ribs_loc[i] >= y_frac:
            y_ribs1 = ribs_loc[i-1]
            #print(y_ribs1)
            y_ribs2 = ribs_loc[i]
            break
    ribs_pitch = (y_ribs2-y_ribs1)*(span/2)
    b = min(ribs_pitch,spar_height)
    a = max(ribs_pitch,spar_height)
    ks_value = np.interp(a/b,ks['ratio'],ks['ks'])
#    print(a/b,ribs_pitch,ks_value,y_ribs1,y_ribs2)
    tau_cr = (np.pi**2*ks_value*E_web)/(12*(1-web_poisson_ratio**2))*(t_spar/b)**2

    return [ribs_pitch,spar_height,a/b,ks_value,tau_cr,shear_strength_web]


#print(get_web_buckling(7,0.4,38))
    
"""
Obtain skin buckling stress
------------------------------------------------------------------------------
"""
kc = pd.read_excel('wing_box_design.xlsx','kc_ss_c')

E_skin = 71*10**9
comp_strength_skin = 441*10**6
skin_poisson_ratio = 0.33

def get_skin_buckling(y_loc,stringer_space, span):
    y_frac = y_loc/(span/2)
    i = 0
    for i in range(1,len(ribs_loc)):
        if ribs_loc[i] >= y_frac:
            y_ribs1 = ribs_loc[i-1]
            #print(y_ribs1)
            y_ribs2 = ribs_loc[i]
            break
    ribs_pitch = (y_ribs2-y_ribs1)*(span/2)
    b = min(ribs_pitch,stringer_space)
    a = max(ribs_pitch,stringer_space)
    kc_value = np.interp(a/b,kc['ratio'],kc['kc'])
    sigma_cr = (np.pi**2*kc_value*E_skin)/(12*(1-skin_poisson_ratio**2))*(t_skin/b)**2
    
    return [ribs_pitch,stringer_space,a/b,kc_value,sigma_cr,comp_strength_skin]

"""
Obtain cross-section properties
y_loc,wing box area,full area,fuel area,moment of inertia,mass at cross section,limit stress 
------------------------------------------------------------------------------
"""

def get_cross_sec_prop(chord,y_loc,span):
    N = 1000
    front_spar = spar_front*chord
    rear_spar = spar_rear*chord
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
    h_front_spar = upper_frontspar_z-lower_frontspar_z
    h_rear_spar = upper_rearspar_z-lower_rearspar_z
    h_ratio = h_front_spar/h_rear_spar

    #-----------------------------Obtain total area and neutral_axis---------------
    tot_stringer_area = upper_stringer*stringer_area_upper+lower_stringer*stringer_area_lower
    front_spar_area = (upper_frontspar_z-lower_frontspar_z)*t_spar
    rear_spar_area = (upper_rearspar_z-lower_rearspar_z)*t_spar
    
    front_spar_area_z = front_spar_area * (upper_frontspar_z+lower_frontspar_z)/2
    rear_spar_area_z = rear_spar_area * (upper_rearspar_z+lower_rearspar_z)/2
    
    skin_box_area = 0
    skin_box_area_z = 0
    skin_box_area_x = 0
    enclosed_area = 0
    for i in range(N):
        x1 = front_spar+i*step_size
        x2 = front_spar+i*step_size+step_size
        z1_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z2_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z1_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        z2_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        
        skin_box_area += abs(x2-x1)*t_skin + abs(x2-x1)*t_skin
        skin_box_area_z += abs(x2-x1)*t_skin*(z2_up+z1_up)/2+\
                            abs(x2-x1)*t_skin*(z2_low+z1_low)/2
        skin_box_area_x += abs(x2-x1)*t_skin*(x2+x1)/2+\
                            abs(x2-x1)*t_skin*(x2+x1)/2
        enclosed_area += abs(((z2_up+z1_up)/2 - (z2_low+z1_low)/2)*(x2-x1))
        
        
    tot_stringer_area_z = 0
    tot_stringer_area_x = 0
    for j in range(upper_stringer):
        xi = stringer_interval_upper*(j+1)+front_spar
        zi = np.interp(xi,foil_upper[:,0],foil_upper[0:,1])
        tot_stringer_area_z += stringer_area_upper*zi
        tot_stringer_area_x += stringer_area_upper*xi
        x_lst.append(xi)
        z_lst.append(zi)
        
        
    for k in range(lower_stringer):
        xi = stringer_interval_lower*(k+1)+front_spar
        zi = np.interp(xi,foil_lower[:,0],foil_lower[:,1])
        tot_stringer_area_z += stringer_area_lower*zi
        tot_stringer_area_z += stringer_area_lower*zi
        x_lst.append(xi)
        z_lst.append(zi)
        
    total_area = skin_box_area+tot_stringer_area+front_spar_area+rear_spar_area
    total_area_z = skin_box_area_z+tot_stringer_area_z+front_spar_area_z+rear_spar_area_z
    
    front_spar_area_x = front_spar_area*front_spar
    rear_spar_area_x = rear_spar_area*rear_spar
    
    total_area_x = skin_box_area_x+tot_stringer_area_x+front_spar_area_x+rear_spar_area_x
    neutral_x = (total_area_x/total_area)
    fuel_area = enclosed_area - total_area
    neutral_z = total_area_z/total_area
    max_z = max(foil_upper[:,1])-neutral_z
    min_z = min(foil_lower[:,1])-neutral_z
    
    #-------------------------Calculate area moment of inertia-----------------
    front_spar_moi = (1/12)*((upper_frontspar_z-lower_frontspar_z)**3)*t_spar + \
                        front_spar_area*((upper_frontspar_z+lower_frontspar_z)/2-
                                         neutral_z)**2
    front_spar_Q = (upper_frontspar_z-lower_frontspar_z)*t_spar*(upper_frontspar_z-
                                         neutral_z)/2
    rear_spar_moi = (1/12)*((upper_rearspar_z-lower_rearspar_z)**3)*t_spar + \
                        rear_spar_area*((upper_rearspar_z+lower_rearspar_z)/2-
                                         neutral_z)**2
    rear_spar_Q = (upper_rearspar_z-lower_rearspar_z)*t_spar*(upper_rearspar_z-
                                         neutral_z)/2

    skin_moi = 0
    upper_skin_Q = 0 #First area moment
    lower_skin_Q = 0
    for i in range(N):
        x1 = front_spar+i*step_size
        x2 = front_spar+i*step_size+step_size
        z1_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z2_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z1_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        z2_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        skin_moi += abs(x2-x1)*t_skin*((z2_up+z1_up)/2-neutral_z)**2
        skin_moi += abs(x2-x1)*t_skin*((z2_low+z1_low)/2-neutral_z)**2
        if x1 < neutral_x:
            upper_skin_Q += abs(x2-x1)*t_skin*((z2_up+z1_up)/2-neutral_z)
        if x1 < neutral_x:
            lower_skin_Q += abs(x2-x1)*t_skin*((z2_low+z1_low)/2-neutral_z)
        
    upper_stringer_moi = 0
    upper_stringer_Q = 0
    for j in range(upper_stringer):
        xi = stringer_interval_upper*(j+1)+front_spar
        zi = np.interp(xi,foil_upper[:,0],foil_upper[:,1])
        upper_stringer_moi += stringer_area_upper*(zi-neutral_z)**2
        if xi < neutral_x:
            upper_stringer_Q += stringer_area_upper*(zi-neutral_z)
    
    lower_stringer_moi = 0
    lower_stringer_Q = 0
    for k in range(lower_stringer):
        xi = stringer_interval_lower*(k+1)+front_spar
        zi = np.interp(xi,foil_lower[:,0],foil_lower[:,1])
        lower_stringer_moi += stringer_area_lower*(zi-neutral_z)**2
        if xi < neutral_x:
            lower_stringer_Q += stringer_area_lower*(zi-neutral_z)

    moi = front_spar_moi+rear_spar_moi+skin_moi+upper_stringer_moi+lower_stringer_moi
    Q = front_spar_Q+upper_stringer_Q+upper_skin_Q
    
    spar_height = h_front_spar
    stringer_space = stringer_interval_upper
    
#    # First moment verification by calculating it at every corner
#    print(Q,rear_spar_Q,lower_stringer_Q,lower_skin_Q)
#    
#    # Verification with CATIA Model
#    print("Verification--------------------------------------------------")
#    print("CATIA Model | Chord = 5m, 2 stringers top and 2 stringers at bottom")
#    print("")
#    print("Front Spar height diff [%]= ", abs(1-(upper_frontspar_z-lower_frontspar_z)/0.688)*100)
#    print("Rear Spar Height diff [%]= ", abs(1-(upper_rearspar_z-lower_rearspar_z)/0.337)*100)
#    print("Cross-section Area diff [%] =",1-(total_area/1.333)*100)
#    print("Neutral-axis [%] =", abs(1-neutral_z/0.1518)*100)
#    print("Moment of inertia [%] =", (1-moi/0.003)*100)
    
    return [y_loc,enclosed_area,total_area,fuel_area,neutral_x,neutral_z,min_z,max_z,moi,\
            spar_height,Q,stringer_space]

##Verification
#stuff = get_cross_sec_prop(5,1,1)
#print(stuff)
    

    

    