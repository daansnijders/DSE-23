# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 22:27:28 2019

@author: Lisa
"""
import numpy as np
import matplotlib.pyplot as plt

from inputs.concept_1 import *
from inputs.constants import *
from inputs.performance_inputs import *


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

    def loading_diagrams(self):
        self.xcg1 = [self.x_cg]
        self.xcg2 = [self.x_cg]
        self.weight = [self.OEW]
          
        self.M_payload_cargo_list = []
        
        for i in range(len(self.x_cargo)):
            self.M_payload_cargo_list.append(self.M_payload_cargo/len(self.x_cargo))
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.x_cargo[i]*self.M_payload_cargo_list[i])/(self.M_payload_cargo_list[i]+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.x_cargo[len(self.x_cargo)-1-i]*self.M_payload_cargo_list[i])/(self.M_payload_cargo_list[i]+self.weight[-1]))
            self.weight.append(self.M_payload_cargo_list[i]+self.weight[-1])
        
        passenger_cg1= np.arange(self.Xfirst, self.Xlast+self.seat_pitch,self.seat_pitch)
        passenger_cg2= np.arange(self.Xlast, self.Xfirst-self.seat_pitch,-self.seat_pitch)
        #print(passenger_cg1)
        
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
        for i in range(len(self.x_fuel)):
            self.xcg1.append((self.weight[-1]*self.xcg1[-1]+self.x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.xcg2.append((self.weight[-1]*self.xcg2[-1]+self.x_fuel[i]*self.M_fuel)/(self.M_fuel+self.weight[-1]))
            self.weight.append(self.M_fuel+self.weight[-1])      
            
            
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

        return min(self.xcg1), max(self.xcg2), max(self.weight)


