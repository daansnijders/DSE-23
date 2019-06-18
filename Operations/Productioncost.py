#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:51:32 2019

@author: coenminderhoud
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt

M_wing= 4000
M_Empennage=200
M_Fuselage=10000
M_Landinggear=500
M_Engines=5000
M_Systems=1500
M_Payloads=9000
M_OEW=35000

WingL=609*2.20462262 * 1.419106771007
WingM=204*2.20462262 * 1.419106771007
WingO=88*2.20462262 * 1.419106771007
EmpennageL=1614*2.20462262 * 1.419106771007
EmpennageM=484*2.20462262 * 1.419106771007
EmpennageO=233*2.20462262 * 1.419106771007
FuselageL=679*2.20462262 * 1.419106771007
FuselageM=484*2.20462262 * 1.419106771007
FuselageO=233*2.20462262 * 1.419106771007
LandinggearL=107*2.20462262 * 1.419106771007
LandinggearM=98*2.20462262 * 1.419106771007
LandinggearO=16*2.20462262 * 1.419106771007
EnginesL=248*2.20462262 * 1.419106771007
EnginesM=91*2.20462262 * 1.419106771007
EnginesO=36*2.20462262 * 1.419106771007
SystemsL=315*2.20462262 * 1.419106771007
SystemsM=91*2.20462262 * 1.419106771007
SystemsO=46*2.20462262 * 1.419106771007
PayloadL=405*2.20462262 * 1.419106771007
PayloadM=100*2.20462262 * 1.419106771007
PayloadO=59*2.20462262 * 1.419106771007
FinalL=58*2.20462262 * 1.419106771007
FinalM=4*2.20462262 * 1.419106771007
FinalO=3*2.20462262 * 1.419106771007

####### calculate cost #######
def findcost(Labor,Materials,Other,Mass):
    C_L=Labor*Mass
    C_M=Materials*Mass
    C_O=Other*Mass
    C_T=C_L+C_M+C_O
    print ('Total cost', C_T ,'Labor cost' ,C_L, 'Material cost' ,C_M, 'Other cost' ,C_O)
    return C_T, C_L, C_M, C_O
###### Wing ######
C_Tw,C_Lw, C_Mw, C_Ow =findcost(WingL,WingM,WingO,M_wing)
###### Empennage ######
C_Te,C_Le, C_Me, C_Oe=findcost(EmpennageL,EmpennageM,EmpennageO,M_Empennage)
###### Fuselage ######
C_Tf,C_Lf, C_Mf, C_Of=findcost(FuselageL,FuselageM,FuselageO,M_Fuselage)
######  Landinggear ######
C_Tl,C_Ll, C_Ml, C_Ol=findcost(LandinggearL,LandinggearM,LandinggearO,M_Landinggear)
###### Engines ######
C_Ten,C_Len, C_Men, C_Oen=findcost(EnginesL,EnginesM,EnginesO,M_Engines)
###### Systems ######
C_Ts,C_Ls, C_Ms, C_Os=findcost(SystemsL,SystemsM,SystemsO,M_Systems)
###### Payloads ######
C_Tp,C_Lp, C_Mp, C_Op=findcost(PayloadL,PayloadM,PayloadO,M_Payloads)
###### Assembly ######
C_TF,C_LF, C_MF, C_OF=findcost(FinalL,FinalM,FinalO,M_OEW)
###### Total ######
totalcost=C_Tw+C_Te+C_Tf+C_Tl+C_Ten+C_Ts+C_Tp+C_TF
print('the total cost is', totalcost)
totalLcost = C_Lw+C_Le+C_Lf+C_Ll+C_Len+C_Ls+C_Lp+C_LF
totalMcost = C_Mw+C_Me+C_Mf+C_Ml+C_Men+C_Ms+C_Mp+C_MF
totalOcost = C_Ow+C_Oe+C_Of+C_Ol+C_Oen+C_Os+C_Op+C_OF



#==============================================================================
#                      Application of learning curve
#==============================================================================
no_of_units = 400
s_labour = 0.85
s_material = 0.95
s_other = 0.95 
n = np.linspace(0, no_of_units, no_of_units)
unit_cost = []
for i in range(no_of_units):
    labour_cost = totalLcost * (i+1)**(np.log(s_labour)/np.log(2.0))
    material_cost = totalMcost * (i+1)**(np.log(s_material)/np.log(2.0))
    other_cost  = totalOcost * (i+1)**(np.log(s_other)/np.log(2.0))
    total_cost = labour_cost + material_cost + other_cost
    unit_cost.append(total_cost)

plt.figure()
plt.plot(n,unit_cost)
plt.show()
















