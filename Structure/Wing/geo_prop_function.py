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

# Spar thickness [m]
t_spar = 0.005

# Skin thickness [m]
t_skin = 0.001

# Front spar location/chord [-]
spar_front = 0.25

# Rear spar location/chord [-]
spar_rear = 0.75

# Chord [m] example
chord = 5


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
    