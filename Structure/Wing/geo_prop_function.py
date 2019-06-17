# -*- coding: utf-8 -*-
"""
Created on Fri May 31 17:27:42 2019
@author: thong

Calculate geometrical property of the cross-section at y_location or chord
of the wing.
"""
import numpy as np
import pandas as pd
import lift_surface_function as lsf
"""
Parameters
------------------------------------------------------------------------------
"""

# Spar thickness [m]
t_spar = 0.015

# Skin thickness [m]
t_skin = 0.003

# Front spar location/chord [-]
spar_front = 0.15#0.1

# Rear spar location/chord [-]
spar_rear = 0.65

# Chord [m] example
chord = 5

stringer_area_upper = 400/10**6
stringer_area_lower = 400/10**6

# Spar materials properties: Aluminium 7150T7751
E_spar = 71*10**9
shear_strength_spar = 462*10**6
spar_poisson_ratio = 0.33
density_spar = 2810
cost_spar = 4.43

# Skin Material: Aluminium 2224A,T351 
E_skin = 72*10**9                 # Pa
comp_strength_skin = 344*10**6      # Pa
skin_poisson_ratio = 0.33
density_skin = 2360#2800                 # kgm^-3
cost_skin = 415#2.8

# Stringer Material: Aluminium 7150T7751
E_stringer = 71*10**9
shear_strength_web = 462*10**6
web_poisson_ratio = 0.33
density_stringer = 2810
cost_stringer = 4.43

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
    tau_cr = (np.pi**2*ks_value*E_spar)/(12*(1-web_poisson_ratio**2))*(t_spar/b)**2

    return [ribs_pitch,spar_height,a/b,ks_value,tau_cr,shear_strength_spar]


#print(get_web_buckling(7,0.4,38))
    
"""
Obtain skin buckling stress
------------------------------------------------------------------------------
"""
kc = pd.read_excel('wing_box_design.xlsx','kc_ss_c')

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
    max_x = rear_spar-neutral_x
    min_x = front_spar-neutral_x
    #-------------------------Average material properties----------------------
    E_modulus = (skin_box_area*E_skin+tot_stringer_area*E_stringer+\
                (front_spar_area+rear_spar_area)*E_spar)/total_area
                
    density = (skin_box_area*density_skin+tot_stringer_area*density_stringer+\
              (front_spar_area+rear_spar_area)*density_spar)/total_area
 
    cost = (skin_box_area*cost_skin+tot_stringer_area*cost_stringer+\
              (front_spar_area+rear_spar_area)*cost_spar)/total_area
    
    #-------------------------Calculate area moment of inertia-----------------
    front_spar_moi_x = (1/12)*((upper_frontspar_z-lower_frontspar_z)**3)*t_spar + \
                        front_spar_area*((upper_frontspar_z+lower_frontspar_z)/2-
                                         neutral_z)**2
    front_spar_moi_z = (1/12)*((upper_frontspar_z-lower_frontspar_z))*t_spar**2 + \
                        front_spar_area*(front_spar-neutral_x)**2

    front_spar_Q_z = (upper_frontspar_z-lower_frontspar_z)*t_spar*(upper_frontspar_z-
                                         neutral_z)/2

    front_spar_Q_x = (upper_frontspar_z-lower_frontspar_z)*t_spar*(front_spar-
                                         neutral_x)

    rear_spar_moi_x = (1/12)*((upper_rearspar_z-lower_rearspar_z)**3)*t_spar + \
                        rear_spar_area*((upper_rearspar_z+lower_rearspar_z)/2-
                                         neutral_z)**2
    rear_spar_moi_z = (1/12)*((upper_rearspar_z-lower_rearspar_z))*t_spar**3 + \
                        rear_spar_area*(rear_spar-
                                         neutral_x)**2

    rear_spar_Q_z = (upper_rearspar_z-lower_rearspar_z)*t_spar*(upper_rearspar_z-
                                         neutral_z)/2

    rear_spar_Q_x = (upper_frontspar_z-lower_frontspar_z)*t_spar*(rear_spar-
                                         neutral_x)
    
    skin_moi_x = 0
    skin_moi_z = 0
    upper_skin_Q_z = 0 #First area moment
    upper_skin_Q_x = 0
    lower_skin_Q_z = 0
    lower_skin_Q_x = 0
    for i in range(N):
        x1 = front_spar+i*step_size
        x2 = front_spar+i*step_size+step_size
        skin_area_i = abs(x2-x1)*t_skin
        z1_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z2_up = np.interp(x1,foil_upper[:,0],foil_upper[:,1])
        z1_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        z2_low = np.interp(x1,foil_lower[:,0],foil_lower[:,1])
        skin_moi_x += skin_area_i*((z2_up+z1_up)/2-neutral_z)**2
        skin_moi_x += skin_area_i*((z2_low+z1_low)/2-neutral_z)**2
        skin_moi_z += skin_area_i*2*(x1-neutral_x)**2
        if x1 < neutral_x:
            upper_skin_Q_z += skin_area_i*((z2_up+z1_up)/2-neutral_z)
            upper_skin_Q_x += skin_area_i*(x1-neutral_x)
        if x1 < neutral_x:
            lower_skin_Q_z += skin_area_i*((z2_low+z1_low)/2-neutral_z)
            lower_skin_Q_x += skin_area_i*(x1-neutral_x)
        
    upper_stringer_moi_x = 0
    upper_stringer_moi_z = 0
    upper_stringer_Q_z = 0
    upper_stringer_Q_x = 0
    for j in range(upper_stringer):
        xi = stringer_interval_upper*(j+1)+front_spar
        zi = np.interp(xi,foil_upper[:,0],foil_upper[:,1])
        upper_stringer_moi_x += stringer_area_upper*(zi-neutral_z)**2
        upper_stringer_moi_z += stringer_area_upper*(xi-neutral_x)**2
        if xi < neutral_x:
            upper_stringer_Q_z += stringer_area_upper*(zi-neutral_z)
            upper_stringer_Q_x += stringer_area_upper*(xi-neutral_x)
    
    lower_stringer_moi_x = 0
    lower_stringer_moi_z = 0
    lower_stringer_Q_z = 0
    lower_stringer_Q_x = 0
    for k in range(lower_stringer):
        xi = stringer_interval_lower*(k+1)+front_spar
        zi = np.interp(xi,foil_lower[:,0],foil_lower[:,1])
        lower_stringer_moi_x += stringer_area_lower*(zi-neutral_z)**2
        lower_stringer_moi_z += stringer_area_lower*(xi-neutral_x)**2
        if xi < neutral_x:
            lower_stringer_Q_z += stringer_area_lower*(zi-neutral_z)
            lower_stringer_Q_x += stringer_area_lower*(xi-neutral_x)
            
    moi_x = front_spar_moi_x+rear_spar_moi_x+skin_moi_x+upper_stringer_moi_x+lower_stringer_moi_x
    moi_z = front_spar_moi_z+rear_spar_moi_z+skin_moi_z+upper_stringer_moi_z+lower_stringer_moi_z

    Q_z = front_spar_Q_z+upper_stringer_Q_z+upper_skin_Q_z
    Q_x = front_spar_Q_x+upper_stringer_Q_x+upper_skin_Q_x
    
    spar_height = h_rear_spar
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
    
    return [y_loc,enclosed_area,total_area,fuel_area,neutral_x,neutral_z,min_z,max_z,moi_x,\
            spar_height,Q_z,stringer_space,moi_z,Q_x,min_x,max_x,E_modulus,density,cost]

##Verification
#stuff = get_cross_sec_prop(5,1,1)
#print(stuff)
    
"""
Get structural weight at every point of the wing
[xi,yi,zi,Fx,Fy,Fz]
"""
def get_struc_force(cross_section_lst,load_factor,sweep_LE,b_wing,Cr,Ct):
    step_size = b_wing/(2*len(cross_section_lst))
    struc_force_lst = []
    for i in range(len(cross_section_lst)):
        yi = cross_section_lst[i,0]
        xi = yi*np.tan(sweep_LE/180*np.pi)+cross_section_lst[i,4]
        zi = 0
        Fx = 0
        Fy = 0
        Fz = -9.81*load_factor*cross_section_lst[i,17]*cross_section_lst[i,2]*step_size
        
        struc_force_lst.append([xi,yi,zi,Fx,Fy,Fz])
    return struc_force_lst
    
"""
Get total cost of the wing
"""
def get_wingbox_mass_cost(cross_section_lst,b_wing):
    step_size = b_wing/(2*len(cross_section_lst))
    mass = 0
    cost = 0
    for i in range(len(cross_section_lst)):
        mass_i = cross_section_lst[i,2]*cross_section_lst[i,17]*step_size
        mass += mass_i
        cost+= mass_i*cross_section_lst[i,18]
    return [mass,cost]
    