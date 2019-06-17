#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:51:32 2019

@author: coenminderhoud
"""
from math import *

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
    C_O=WingO*Mass
    C_T=C_L+C_M+C_O
    print ('Total cost', C_T ,'Labor cost' ,C_L, 'Material cost' ,C_M, 'Other cost' ,C_O)
    return C_T
###### Wing ######
C_Tw=findcost(WingL,WingM,WingO,M_wing)
###### Empennage ######
C_Te=findcost(EmpennageL,EmpennageM,EmpennageO,M_Empennage)
###### Fuselage ######
C_Tf=findcost(FuselageL,FuselageM,FuselageO,M_Fuselage)
######  Landinggear ######
C_Tl=findcost(LandinggearL,LandinggearM,LandinggearO,M_Landinggear)
###### Engines ######
C_Ten=findcost(EnginesL,EnginesM,EnginesO,M_Engines)
###### Systems ######
C_Ts=findcost(SystemsL,SystemsM,SystemsO,M_Systems)
###### Payloads ######
C_Tp=findcost(PayloadL,PayloadM,PayloadO,M_Payloads)
###### Assembly ######
C_Tf=findcost(FinalL,FinalM,FinalO,M_OEW)
###### Total ######
totalcost=C_Tw+C_Te+C_Tf+C_Tl+C_Tl+C_Ten+C_Ts+C_Tp+C_Tf
print('the total cost is', totalcost)