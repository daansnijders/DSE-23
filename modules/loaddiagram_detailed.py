# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 22:27:28 2019

@author: Lisa
"""
import numpy as np
import matplotlib.pyplot as plt

from inputs.concept_1 import *
from inputs.constants import *


"Values to input once known"
x_cg_wing = 17.
x_cargo = [5, 20]
x_fuel = [x_cg_wing]
Xfirst = 5

class Loading_diagram:
    def __init__(self, l_f, l_cabin, seat_pitch, N_pax, N_sa, OEW, MTOW, M_MZF,x_cg, y_cg, z_cg, MAC, S, b, A, Xfirst, M_payload, M_cargo_available, M_fuel, M_pax, M_carry_on):
        self.l_f=l_f
        self.l_cabin=l_cabin
        self.seat_pitch=seat_pitch
        self.N_pax=N_pax
        self.N_sa=N_sa
        self.OEW=OEW
        self.MTOW=MTOW
        self.M_MZF=M_MZF
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
        self.M_pax=M_pax
        self.M_carry_on=M_carry_on

        self.N_rows = self.N_pax/self.N_sa
        self.Xlast=self.Xfirst+seat_pitch*(self.N_rows-1)
        self.M_pax_cabin = M_pax + M_carry_on
        self.M_pass_total=self.M_pax_cabin*self.N_pax
    
        self.M_payload_cargo=self.M_payload-self.M_pass_total


    def loading_diagrams(self):
        self.xcg1 = [self.x_cg]
        self.xcg2 = [self.x_cg]
        self.weight = [self.OEW]
          
        self.M_payload_cargo_list = []
        self.M_fuel_list = []
        
        
        for i in range(len(x_cargo)):
            self.M_payload_cargo_list.append(self.M_payload_cargo/len(x_cargo))
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+x_cargo[i]*self.M_payload_cargo_list[i])/(self.M_payload_cargo_list[i]+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+x_cargo[len(x_cargo)-1-i]*self.M_payload_cargo_list[i])/(self.M_payload_cargo_list[i]+self.weight[-1]))
            self.weight.append(self.M_payload_cargo_list[i]+self.weight[-1])
        
        passenger_cg1= np.arange(self.Xfirst, self.Xlast+self.seat_pitch,self.seat_pitch)
        passenger_cg2= np.arange(self.Xlast, self.Xfirst-self.seat_pitch,-self.seat_pitch)
        
        for i in range(len(passenger_cg1)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+passenger_cg1[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+passenger_cg2[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))        
            self.weight.append(2*self.M_pax_cabin+self.weight[-1])

        for i in range(len(passenger_cg1)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+passenger_cg1[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+passenger_cg2[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))        
            self.weight.append(2*self.M_pax_cabin+self.weight[-1])
            
        for i in range(len(passenger_cg1)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+passenger_cg1[i]*self.M_pax_cabin)/(self.M_pax_cabin+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+passenger_cg2[i]*self.M_pax_cabin)/(self.M_pax_cabin+self.weight[-1]))        
            self.weight.append(self.M_pax_cabin+self.weight[-1])            
        
        "Only considered one fuel tank for now at place of wing"        
        for i in range(len(x_fuel)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.weight.append(self.M_fuel+self.weight[-1])      
        
        print(min(self.xcg1))           
        
        plt.figure()   
        plt.plot(self.xcg1, self.weight, color='blue', marker='o')
        plt.plot(self.xcg2, self.weight, color='green', marker='o')
        #plt.hlines(Weight[23],min(xcg), max(xcg),'r')
        #plt.hlines(Weight[45],min(xcg), max(xcg), 'r')
        #plt.hlines(self.weight[-1],min(self.xcg1), max(self.xcg2), 'r')
        #plt.hlines(self.weight[-2],min(self.xcg1), max(self.xcg2), 'r')
        plt.vlines(min(self.xcg1)-0.02,self.weight[0], self.weight[-1], 'k')
        plt.vlines(max(self.xcg2)+0.02,self.weight[0], self.weight[-1], 'k')
        plt.title('Title', fontsize=14)
        plt.xlabel('Center of gravity location from nose [m]', fontsize=12)
        plt.ylabel('Mass[kg]', fontsize=12)
        plt.show()
        
        
        return self.xcg1[-1], self.xcg2[-1], self.weight[-1]
##
##    
##    #Xlemac_plus=(Xlemac+0.1*fuselage_length)
##    #Xlemac_min=Xlemac-0.1*fuselage_length
##    #
##    #xcg1, xcg2,Weight1= loading_diagrams(Xlemac)
##    #xcg1_plus, xcg2_plus,Weight2= loading_diagrams(Xlemac_plus)
##    #xcg1_min, xcg2_min,Weigth3= loading_diagrams(Xlemac_min)
##    
##    
##    #min_cg=[min(xcg1_min)-2,min(xcg1)-2,min(xcg1_plus)-2]    
##    #max_cg=[max(xcg2_min)+2,max(xcg2)+2,max(xcg2_plus)+2]        
##    #long_pos=[ Xlemac_min/fuselage_length,Xlemac/fuselage_length, Xlemac_plus/fuselage_length]
##    #wing longitudinal position vs cg position 
##    
##    
##    #plt.figure()
##    #plt.plot(min_cg,long_pos,label='Most forward CG location')
##    #plt.plot(max_cg,long_pos,label='Most aft CG location')   
##    #plt.xlabel('center of gravity range [% of MAC]', fontsize=12) 
##    #plt.ylabel('Longitudinal wing position (Xlemac/Lfuselage) [-]', fontsize=12)
##    #plt.title('c.g. range vs. wing longitudinal position diagram ',fontsize=14)
##    #plt.legend()
##    ##plt.show()
##
##
##
##
##
##
##
##
##
##
##
