# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:17:48 2019

@author: daansnijders
"""

from inputs.concept_1 import *
from inputs.constants import *


class check_ground():
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
            
        assert min(self.frac) > 0.08
        assert max(self.frac) < 0.20
            
        return min(self.frac), max(self.frac), self.frac
    
    
    
def update_x_mlg(z_cg,theta_rad, beta_rad, x_cg, stroke, l_f):
    beta_rad_correct = beta_rad- np.pi/2
    
    scrape = tangent(l_f,z_cg,theta_rad)
    tip_over = tangent(x_cg,z_cg,beta_rad_correct)
    
    x_mlg = (tip_over[1] - scrape[1] + stroke)/(scrape[0] - (tip_over[0])) 
    
    return x_mlg

def update_z_mlg(x_mlg,beta_rad,x_cg, z_cg):
    beta_rad_correct = beta_rad - np.pi/2

    tip_over = tangent(x_cg,z_cg,beta_rad_correct)
    
    
    z_mlg = tip_over[0] * x_mlg + tip_over[1]
    return z_mlg

def update_l_mw(x_mlg,x_cg):
    return x_mlg - x_cg

def update_l_nw(l_w,P_mw,N_mw,P_nw,N_nw):
    return [(l_w[i]*P_mw[i]*N_mw)/(P_nw[i]*N_nw) for i in range(3)]

def update_y_mlg(z_cg,z_mlg,l_n,l_w):
    z_t = b/(2) * np.tan(dihedral_rad)
    
    y_mlg1 = [(l_n[i] + l_w[i])/(np.sqrt((l_n[i]**2 \
              *np.tan(psi_rad)**2)/(z_cg - z_mlg)**2 - 1)) for i in range(3)]
    y_mlg2 = [y_engine - (z_engine - z_mlg)/np.tan(phi_rad) for i in range(3)]                           
    y_mlg3 = [b/2 - (z_t-z_mlg)/np.tan(phi_rad) for i in range(3)] 
    
    y_mlg = max([max([y_mlg1[i],y_mlg2[i],y_mlg3[i]]) for i in range(3)])
    assert y_mlg < 9/2
    return y_mlg
    
    
    