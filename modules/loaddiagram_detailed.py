# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 22:27:28 2019

@author: Lisa
"""
import numpy as np
import matplotlib.pyplot as plt

from inputs.concept1 import *
from inputs.constants import *

x_cargo_front = 8
x_cargo_aft = 20

class Loading_diagram:
    def __init__(self, l_f, l_cabin, seat_pitch, N_pax, N_sa, OEW, MTOW, M_MZF, M_TO, N_rows, x_cg, y_cg, z_cg, MAC, S, b, A, Xfirst, M_payload, M_cargo_available, M_fuel, W_pax, W_carry_on):
        self.l_f=l_f
        self.l_cabin=l_cabin
        self.seat_pitch=seat_pitch
        self.N_pax=N_pax
        self.N_sa=N_sa
        self.OEW=OEW
        self.MTOW=MTOW
        self.M_MZF=M_MZF
        self.M_TO=M_TO
        self.N_rows=N_rows
        self.x_cg=x_cg
        self.y_cg=y_cg
        self.z_cg=z_cg
        self.MAC=MAC
        self.S=S
        self.b=b
        self.A=A
        self.Xfirst=Xfirst
        self.M_payload=M_payload
        self.M_cargo_available=M_cargo_available
        self.M_fuel=M_fuel
        self.W_pax=W_pax
        self.W_carry_on=W_carry_on


    self.Xlast=Xfirst+seat_pitch*(N_rows-1)
    
    self.M_pax_cabin = W_pax + W_carry_on
    self.W_pass_total=W_pax_cabin*N_pax
    
    self.W_payload=self.M_payload_total-self.W_pass_total


    def loading_diagrams(self):
        xcg = [x_cg]
        weight = [OEW]
    
        W_payload_front=W_payload/2
        W_payload_aft=W_payload/2
    
        passenger_cg1= np.arange(Xfirst, Xlast-seatpitch,seatpitch)
        passenger_cg2= np.arange(Xlast, Xfirst-seatpitch,-seatpitch)
        
        print(passenger_cg1,passenger_cg2)
        
        for i in range(len(passenger_cg1)):
            cg_prev1=xcg1[i+2]
            cg_prev2=xcg2[i+2]
            W_old=Weight[i+2]
            W_new=W_old+W_passenger*2
            cg2=(W_passenger*2*passenger_cg2[i]+W_old*cg_prev2)/(W_new)
            cg1=(W_passenger*2*passenger_cg1[i]+W_old*cg_prev1)/(W_new)
            xcg1.append(cg1)
            xcg2.append(cg2)
            Weight.append(W_new)
            
            
         
        for n in range(len(passenger_cg1), (len(passenger_cg1))*2):
            cg_prev1=xcg1[n+2]
            cg_prev2=xcg2[n+2]
            W_old=Weight[n+2]
            W_new=W_old+W_passenger*2
            cg2=(W_passenger*2*passenger_cg2[n-22]+W_old*cg_prev2)/(W_new)
            cg1=(W_passenger*2*passenger_cg1[n-22]+W_old*cg_prev1)/(W_new)
            xcg1.append(cg1)
            xcg2.append(cg2)
            Weight.append(W_new)   
            
        for p in range(len(passenger_cg1)*2, (len(passenger_cg1))*3-1):
            cg_prev1=xcg1[p+2]
            cg_prev2=xcg2[p+2]
            W_old=Weight[p+2]
            W_new=W_old+W_passenger
            cg2=(W_passenger*passenger_cg2[p+1-44]+W_old*cg_prev2)/(W_new)
            cg1=(W_passenger*passenger_cg1[p-44]+W_old*cg_prev1)/(W_new)
            xcg1.append(cg1)
            xcg2.append(cg2)
            Weight.append(W_new)   
        
        
        Wtot=W_fuel+Weight[-1]   
        
        cgfuel=(Weight[-1]*xcg1[-1]+W_fuel*cg_fuel)/Wtot 
        xcg1.append(cgfuel)
        xcg2.append(cgfuel)
        Weight.append(Wtot)
        
        xcg1=np.array(xcg1)*100
        xcg2=np.array(xcg2)*100
        
        plt.figure()   
        plt.plot(xcg1, Weight, color='blue', marker='o')
        plt.plot(xcg2, Weight, color='green', marker='o')
        #plt.hlines(Weight[23],min(xcg), max(xcg),'r')
        #plt.hlines(Weight[45],min(xcg), max(xcg), 'r')
        plt.hlines(Weight[-1],min(xcg1), max(xcg2), 'r')
        plt.hlines(Weight[-2],min(xcg1), max(xcg2), 'r')
        plt.vlines(min(xcg1)-2,Weight[0], Weight[-1], 'k')
        plt.vlines(max(xcg2)+2,Weight[0], Weight[-1], 'k')
        plt.title('Load diagram for Xlemac at %2.3f m ' %(Xlemac), fontsize=14)
        plt.xlabel('Center of gravity location[% of MAC]', fontsize=12)
        plt.ylabel('Mass[kg]', fontsize=12)
        plt.show()
        
        return xcg1, xcg2,Weight

    
    #Xlemac_plus=(Xlemac+0.1*fuselage_length)
    #Xlemac_min=Xlemac-0.1*fuselage_length
    #
    #xcg1, xcg2,Weight1= loading_diagrams(Xlemac)
    #xcg1_plus, xcg2_plus,Weight2= loading_diagrams(Xlemac_plus)
    #xcg1_min, xcg2_min,Weigth3= loading_diagrams(Xlemac_min)
    
    
    #min_cg=[min(xcg1_min)-2,min(xcg1)-2,min(xcg1_plus)-2]    
    #max_cg=[max(xcg2_min)+2,max(xcg2)+2,max(xcg2_plus)+2]        
    #long_pos=[ Xlemac_min/fuselage_length,Xlemac/fuselage_length, Xlemac_plus/fuselage_length]
    #wing longitudinal position vs cg position 
    
    
    #plt.figure()
    #plt.plot(min_cg,long_pos,label='Most forward CG location')
    #plt.plot(max_cg,long_pos,label='Most aft CG location')   
    #plt.xlabel('center of gravity range [% of MAC]', fontsize=12) 
    #plt.ylabel('Longitudinal wing position (Xlemac/Lfuselage) [-]', fontsize=12)
    #plt.title('c.g. range vs. wing longitudinal position diagram ',fontsize=14)
    #plt.legend()
    ##plt.show()











