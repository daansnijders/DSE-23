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


class Loading_diagram:
    def __init__(self, x_cargo, l_f, l_cabin, seat_pitch, N_pax, N_sa, OEW, MTOW, x_cg, y_cg, z_cg, MAC, S, b, A, Xfirst, M_payload, M_cargo_available, M_fuel, M_pax, M_carry_on, x_cg_wing, config):
        self.x_cargo=x_cargo
        self.l_f=l_f
        self.l_cabin=l_cabin
        self.seat_pitch=seat_pitch
        self.N_pax=N_pax
        self.N_sa=N_sa
        self.OEW=OEW
        self.MTOW=MTOW
        self.x_cg=x_cg
        self.y_cg=y_cg
        self.x_cg_wing=x_cg_wing
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
        self.M_payload2=N_pax*100.5
        self.N_rows = self.N_pax/self.N_sa
        self.Xlast=self.Xfirst+seat_pitch*(self.N_rows-1)
        self.M_pax_cabin = M_pax + M_carry_on
        self.M_pass_total=self.M_pax_cabin*self.N_pax
        self.M_payload_cargo=self.M_payload2-self.M_pass_total
        self.x_fuel = [x_cg_wing]
        self.config=config

    def loading_diagrams_pass(self):
        self.xcg1 = [self.x_cg]
        self.xcg2 = [self.x_cg]
        self.weight = [self.OEW]
          
        self.M_payload_cargo_list = []
        
        for i in range(len(self.x_cargo)):
            self.M_payload_cargo_list.append(self.M_payload_cargo/len(self.x_cargo))
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.x_cargo[i]*self.M_payload_cargo_list[i])/(self.M_payload_cargo_list[i]+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.x_cargo[len(self.x_cargo)-1-i]*self.M_payload_cargo_list[i])/(self.M_payload_cargo_list[i]+self.weight[-1]))
            self.weight.append(self.M_payload_cargo_list[i]+self.weight[-1])
        
        self.passenger_cg1= np.arange(self.Xfirst, self.Xlast+self.seat_pitch,self.seat_pitch)
        self.passenger_cg2= np.arange(self.Xlast, self.Xfirst-self.seat_pitch,-self.seat_pitch)
        #print(self.passenger_cg1)
        
        for i in range(len(self.passenger_cg1)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.passenger_cg1[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.passenger_cg2[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))        
            self.weight.append(2*self.M_pax_cabin+self.weight[-1])

        for i in range(len(self.passenger_cg1)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.passenger_cg1[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.passenger_cg2[i]*2*self.M_pax_cabin)/(2*self.M_pax_cabin+self.weight[-1]))        
            self.weight.append(2*self.M_pax_cabin+self.weight[-1])
            
        for i in range(len(self.passenger_cg1)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.passenger_cg1[i]*self.M_pax_cabin)/(self.M_pax_cabin+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.passenger_cg2[i]*self.M_pax_cabin)/(self.M_pax_cabin+self.weight[-1]))        
            self.weight.append(self.M_pax_cabin+self.weight[-1])            
        

        "Only considered one fuel tank for now at place of wing"        
        for i in range(len(self.x_fuel)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.weight.append(self.M_fuel+self.weight[-1])      
        

        self.xcg_max = self.xcg1[-1]
        if (self.xcg1[-2]>=self.xcg_max):
            self.xcg_max = self.xcg1[-2]

        self.xcg_min = self.xcg1[-1]
        if (self.xcg1[-2]<=self.xcg_min):
           self.xcg_min = self.xcg1[-2]

        if (self.config == 2 or self.config == 3):
            self.xcg_max = self.xcg_max - (l_cabin[1]-l_cabin[0])
            self.xcg_min = self.xcg_min - (l_cabin[1]-l_cabin[0])
        
        plt.figure()   
        plt.plot(self.xcg1, self.weight, color='blue', marker='o')
        plt.plot(self.xcg2, self.weight, color='green', marker='o')
        #plt.hlines(Weight[23],min(xcg), max(xcg),'r')
        #plt.hlines(Weight[45],min(xcg), max(xcg), 'r')
        #plt.hlines(self.weight[-1],min(self.xcg1), max(self.xcg2), 'r')
        #plt.hlines(self.weight[-2],min(self.xcg1), max(self.xcg2), 'r')
        plt.vlines(min(self.xcg1)-0.02,self.weight[0], self.weight[-1], 'k')
        plt.vlines(max(self.xcg2)+0.02,self.weight[0], self.weight[-1], 'k')
        plt.title('C.g. location for configuration: %i' %self.config, fontsize=14)
        plt.xlabel('Center of gravity location from nose [m]', fontsize=12)
        plt.ylabel('Mass[kg]', fontsize=12)
        plt.show()

        return self.xcg1, self.xcg2, self.weight, self.xcg_max, self.xcg_min
    
    def loading_diagrams_fuel(self):
        self.xcg1 = [self.x_cg]
        self.xcg2 = [self.x_cg]
        self.weight = [self.OEW]        

        "Only considered one fuel tank for now at place of wing"        
        for i in range(len(self.x_fuel)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.weight.append(self.M_fuel+self.weight[-1])      
            
        self.xcg_max = self.xcg1[-1]
        if (self.xcg1[-2]>=self.xcg_max):
           self.xcg_max = self.xcg1[-2]
        
        self.xcg_min = self.xcg1[-1]
        if (self.xcg1[-2]<=self.xcg_min):
           self.xcg_min = self.xcg1[-2]
        
        if (self.config == 2 or self.config == 3):
            self.xcg_max = self.xcg_max - (l_cabin[1]-l_cabin[0])
            self.xcg_min = self.xcg_min - (l_cabin[1]-l_cabin[0])

        
        plt.figure()   
        plt.plot(self.xcg1, self.weight, color='blue', marker='o')
        plt.plot(self.xcg2, self.weight, color='green', marker='o')
        #plt.hlines(Weight[23],min(xcg), max(xcg),'r')
        #plt.hlines(Weight[45],min(xcg), max(xcg), 'r')
        #plt.hlines(self.weight[-1],min(self.xcg1), max(self.xcg2), 'r')
        #plt.hlines(self.weight[-2],min(self.xcg1), max(self.xcg2), 'r')
        plt.vlines(min(self.xcg1)-0.02,self.weight[0], self.weight[-1], 'k')
        plt.vlines(max(self.xcg2)+0.02,self.weight[0], self.weight[-1], 'k')
        plt.title('C.g. location for configuration: %i' %self.config, fontsize=14)
        plt.xlabel('Center of gravity location from nose [m]', fontsize=12)
        plt.ylabel('Mass[kg]', fontsize=12)
        plt.show()

        return self.xcg1, self.xcg2, self.weight, self.xcg_max, self.xcg_min
