# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:17:48 2019

@author: daansnijders
"""

x_nlg = [1., 1., 1.]
x_mlg = [13., 19., 19.]


from inputs.concept_1 import *
from inputs.constants import *


class Stability_check_ground:
    def __init__ (self, cg1_pass, cg2_pass, weight_pass, cg1_fuel, cg2_fuel, weight_fuel, x_nlg, x_mlg):
        self.cg1_pass=cg1_pass
        self.cg2_pass=cg2_pass
        self.weight_pass=weight_pass
        self.cg1_fuel=cg1_fuel
        self.cg2_fuel=cg2_fuel
        self.weight_fuel=weight_fuel
        self.x_nlg=x_nlg
        self.x_mlg=x_mlg
        
        
                
               
    def check_equilibrium(self):
        self.b_n = [0]*(2*len(self.weight_pass))
        self.b_m = [0]*(2*len(self.weight_pass))
        self.F_n = [0]*(2*len(self.weight_pass))
        self.F_m = [0]*(2*len(self.weight_pass))
        self.frac = [0]*(2*len(self.weight_pass))
        
        "Complete loading diagram"
        for i in range(len(self.weight_pass)):
            self.b_n[i] = self.cg1_pass[i] - self.x_nlg
            self.b_m[i] = self.x_mlg - self.cg1_pass[i] 
            self.b_n[(i+len(self.weight_pass))] = self.cg2_pass[i] - self.x_nlg
            self.b_m[(i+len(self.weight_pass))] =  self.x_mlg - self.cg2_pass[i]
                
        
        for i in range(len(self.weight_pass)):
            self.F_n[i] = (self.weight_pass[i]*self.b_m[i])/(self.b_m[i]+self.b_n[i])
            self.F_n[(i+len(self.weight_pass))] = (self.weight_pass[i]*self.b_m[(i+len(self.weight_pass))])/(self.b_m[i]+self.b_n[i])
        
        for i in range(len(self.weight_pass)):
            self.F_m[i] = (self.weight_pass[i]-self.F_n[i])/2
            self.F_m[(i+len(self.weight_pass))] = (self.weight_pass[i]-self.F_n[(i+len(self.weight_pass))])/2
        
        
        for i in range(len(self.weight_pass)):
            self.frac[i] = self.F_n[i]/self.weight_pass[i]
            self.frac[(i+len(self.weight_pass))] = self.F_n[(i+len(self.weight_pass))]/self.weight_pass[i]
        
        "Only OEW with fuel"
        for i in range(len(self.weight_fuel)):
            self.b_n.append(self.cg1_fuel[i] - self.x_nlg)
            self.b_m.append(self.x_mlg - self.cg1_fuel[i])
        for i in range(len(self.weight_fuel)):
            self.b_n.append(self.cg2_fuel[i] - self.x_nlg)
            self.b_m.append(self.x_mlg - self.cg2_fuel[i])
                
        for i in range(len(self.weight_fuel)):
            self.F_n.append((self.weight_fuel[i]*self.b_m[(2*len(self.weight_pass)+i)])/(self.b_m[(2*len(self.weight_pass)+i)]+self.b_n[(2*len(self.weight_pass)+i)]))
        for i in range(len(self.weight_fuel)):
            self.F_n.append((self.weight_fuel[i]*self.b_m[(i+len(self.weight_fuel)+2*len(self.weight_pass))])/(self.b_m[(i+len(self.weight_fuel)+2*len(self.weight_pass))]+self.b_n[(i+len(self.weight_fuel)+2*len(self.weight_pass))]))
        
        for i in range(len(self.weight_fuel)):
            self.F_m.append((self.weight_fuel[i]-self.F_n[(2*len(self.weight_pass)+i)])/2)
        for i in range(len(self.weight_fuel)):
            self.F_m.append((self.weight_fuel[i]-self.F_n[(i+len(self.weight_fuel)+2*len(self.weight_pass))])/2)
        
        for i in range(len(self.weight_fuel)):
            self.frac.append(self.F_n[(2*len(self.weight_pass)+i)]/self.weight_fuel[i])
        for i in range(len(self.weight_fuel)):
            self.frac.append(self.F_n[(i+len(self.weight_fuel)+2*len(self.weight_pass))]/self.weight_fuel[i])
        
                
        return min(self.frac), max(self.frac), self.frac