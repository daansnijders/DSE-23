# -*- coding: utf-8 -*-
"""
Created on Mon May 13 17:29:43 2019

@author: Lisa
"""
from math import *
import numpy as np

#inputs needed

H=37000
M_cr=0.75
M_x=0.935
S=100
MTOW=40000*9.81
A=11


def planformsizing(H,M_cr,S,MTOW,A,):

    M_dd=M_cr+0.03
    if M_cr>0.7:
        sweep_c_4=arccos(0.75*M_x/M_dd)
    else:
        sweep_c_4=1
        
    b=sqrt(S*A)
    taper=0.2*(2-sweep_c_4*pi/180)
    
    c_r=2*S/((1+taper)*b)
    c_t=c_r*taper
    
    dihedral=deg2rad(3-rad2deg(sweep_c_4)/10) #add 2 for low wing , -2 for high wing
    
    p=p0*(1-0.0065*H/(288.15))**(9.81/(287*0.0065))
    q=0.5*1.4*p*M_cr**2
    
    Cl=MTOW/(q*S)
    MAC=c_r*2/3*(1+taper+taper**2)/(1+taper)
    Y_mac=b/2*(c_r-MAC)/(c_r-c_t)
    
    sweep_c_2=arctan(tan(sweep_c_4)-(1-taper/(a*(1+taper))))
    
    t_c=(cos(sweep_c_2)**3*(M_x-M_dd*cos(sweep_c_2))-0.115*Cl**1.5)/cos(sweep_c_2)
    if t_c>0.18:
        t_c=0.18
        
    return b,taper,c_r,c_t,dihedral,Cl, MAC, Y_mac, sweep_c_4, sweep_c_2, t_c

    
