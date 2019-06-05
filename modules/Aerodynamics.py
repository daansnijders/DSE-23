# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np

class HLD_class:
    def __init__(self,Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer):
        self.Cl_land        = Cl_land
        self.Cl_clean       = Cl_clean
        self.S              = S
        self.A              = A
        self.lambda_4_rad   = lambda_4_rad
        self.taper_ratio    = taper_ratio
        self.CL_alpha_clean = CL_alpha
        self.lambda_le_rad  = lambda_le_rad
        self.b              = sqrt(A*S)
        self.Cr             = Cr
        self.d_f_outer      = d_f_outer
    
    def HLD(self):
        Delta_CLmax = self.Cl_land - self.Cl_clean    #[-i + self.Cl_land for i in self.Cl_clean]
        hl = 0.65              # Location hinge line on chord
        lambda_hl_rad = np.arctan(np.tan(self.lambda_4_rad)-(4/self.A)*(hl-1/4)*(1-self.taper_ratio)/(1+self.taper_ratio))
        c_prime = 1 + 0.57*(1-hl)
        Delta_cl = 1.6*(c_prime)
        SWF = (Delta_CLmax *self.S)/(0.9*Delta_cl*np.cos(lambda_hl_rad))
        #[(x *self.S)/(0.9*Delta_cl*np.cos(lambda_hl_rad)) for x in Delta_CLmax]
        SWF_S = SWF/self.S
        delta_alpha = np.deg2rad(-15) * SWF_S * np.cos(lambda_hl_rad)
        
        Sprime_S = 1 + SWF_S * (c_prime - 1)
        
        CL_alpha_flapped = Sprime_S * self.CL_alpha_clean
        
        HLD_clearance = 0.5     #Clearance between fuselage and the HLD's 
        
        """ Calculate span of the flap """
        #x, h1, h2, and h3 are only used for calculation purposes
        x = self.taper_ratio * self.b/2 / (1 - self.taper_ratio)
        h2 = self.b / 2 + x
        S_wet = 0
        h1 = h2 - self.d_f_outer/2 - HLD_clearance
        c_flap_start = h1/h2*self.Cr
        i = 1
        while S_wet <= SWF :
            h3 = h1 - i*0.001
            c_flap_end = h3/h2*self.Cr
            S_wet = 2*((c_flap_start + c_flap_end) / 2 * (i*0.001))
            i += 1
        b_flap = h1 - h3
        
        SWF_LE = (0.1*self.S)/(0.9*0.3*self.lambda_le_rad)
        """ Calculate span of the slat """
        #x, h1, h2, and h3 are only used for calculation purposes
        x = self.taper_ratio * self.b/2 / (1 - self.taper_ratio)
        h2 = self.b / 2 + x
        S_wet = 0
        h1 = h2 - self.d_f_outer/2 - HLD_clearance
        c_slat_start = h1/h2*self.Cr
        i = 1
        while S_wet <= SWF_LE :
            h3 = h1 - i*0.001
            c_slat_end = h3/h2*self.Cr
            S_wet = 2*((c_slat_start + c_slat_end) / 2 * (i*0.001))
            i += 1
        b_slat = h1 - h3
        
        return(SWF, b_flap, SWF_LE, b_slat)
        
class Drag:
    def __init__(self,S,A,rho,rho_0,l_f,V_cruise,V_TO,mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes,CL_alpha,l_cockpit, l_cabin, l_tail,lambda_2_rad):
        self.S              = S
        self.A              = A
        self.rho            = rho
        self.rho_0          = rho_0
        self.l_f            = l_f
        self.V_cruise       = V_cruise
        self.mu_37          = mu_37
        self.MAC            = MAC
        self.Ct             = Ct
        self.Cr             = Cr
        self.b              = b
        self.taper_ratio    = taper_ratio
        self.d_f_outer      = d_f_outer
        self.lambda_le_rad  = lambda_le_rad
        self.V_TO           = V_TO
        self.mu_sl          = mu_sl
        self.M              = 0.75
        self.CLdes          = CLdes
        self.CL_alpha       = CL_alpha
        self.l_cockpit      = l_cockpit
        self.l_cabin        = l_cabin
        self.l_tail         = l_tail
        self.lambda_2_rad   = lambda_2_rad
        
    def wing_drag(self):
        Re_f  = self.rho   * self.V_cruise * self.l_f / self.mu_37
        Re_f0 = self.rho_0 * self.V_TO     * self.l_f / self.mu_sl
        #These result in
        R_wf = 1.01         #(Figure 4.1) 
        
        cos_lambda_2_rad   = cos(lambda_2_rad)
        #This results in
        R_LS   = 1.21        #(Figure 4.2)
        R_LS_c = 1.21        #(Figure 4.2)
        
        L_prime = 2.0
        C_f_w = .00265
        t_c = 0.15
        
        """ Calculate chord at fuselage-wing intersection """
        x = self.taper_ratio * self.b/2 / (1 - self.taper_ratio)
        h2 = self.b / 2 + x
        h1 = h2 - self.d_f_outer/2
        c_fuselage_wing = h1/h2*self.Cr
        
        S_wet = 2*(2*((Ct + c_fuselage_wing) / 2 * (self.b/2 - self.d_f_outer/2)))
        
        C_D_0_W = R_wf * R_LS * C_f_w * (1 + L_prime * (t_c) + 100 * (t_c)**4) * S_wet/self.S
        
        """ C_D_L_w """
        r_LE = 0.687                  #Leading Edge radius
        RE_LER = self.rho * self.V_cruise * r_LE / self.mu_37
        R_par = RE_LER * 1/(tan(self.lambda_le_rad)) * sqrt(1 - (self.M*cos(self.lambda_le_rad))**2)
        R_par2 = self.A * self.taper_ratio / cos(self.lambda_le_rad) 
        #This results in 
        R = 0.95
        
        beta = sqrt(1-self.M**2)
        c_l_alpha = np.deg2rad((1.3 + 0.5)/(7+9))
        k = c_l_alpha/(2*pi / beta)
        C_L_a_w = (2*pi*self.A)/(2 + ((self.A * beta / k)**2 * (1 + (tan(self.lambda_2_rad)))))
        
        e = 1.1*(C_L_a_w / self.A)*(R * (C_L_a_w / self.A) + (1-R)*pi)
        
        C_L_w = 1.05 * self.CLdes
        C_D_L_w = C_L_w**2 / (pi * self.A * e)
        
        C_D_w = C_D_0_w + C_D_L_w
        
        return (C_D_w)
        
                        
    def fuse_drag(self):
        Re_f  = self.rho   * self.V_cruise * self.l_f / self.mu_37
        Re_f0 = self.rho_0 * self.V_TO     * self.l_f / self.mu_sl
        #RE at service ceiling results in (same for all configurations)
        Rwf = 1.01          #Figure 4.1
        Cf_fus = 0.0019     #Figure 4.3
        ratio = self.l_f/self.d_f_outer
        Swet_fus = pi*self.d_f_outer*self.l_f*(1-2/ratio)**(2/3)*(1+1/(ratio)**2)
        
        CD0_fus = Rwf*Cf_fus*(1+60/(self.l_f/self.d_f_outer)**3 + 0.0025*(self.l_f/self.d_f_outer))*Swet_fus/self.S
        
        CL0 = self.CL_alpha*(pi/180)*5.5
        alpha = (self.CLdes -CL0)/self.CL_alpha
        
        eta1 = 0.66         #Figure 4.19
        eta2 = 0.68         #Figure 4.19
        cdc = 1.2           #Figure 4.20
        Splf = 0.5*self.d_f_outer*(self.l_tail+self.l_cockpit) + self.d_f_outer*self.l_cabin
        
        if self.l_cabin == 19.44:
            CDL_fus = eta1*cdc*alpha**3*(Splf/self.S)
        else:
            CDL_fus = eta2*cdc*alpha**3*(Splf/self.S)
        
        CD_fus_sub = CD0_fus + CDL_fus
        
        
        CDf_fus = Cf_fus*(Swet_fus/self.S)
        CDp_fus = Cf_fus*(60/(self.l_f/self.d_f_outer)**3 + 0.0025*(self.l_f/self.d_f_outer))*Swet_fus/self.S
        CD_wave = 0.005     #Figure 4.22
        
        CD_fus_trans = Rwf*(CDf_fus + CDp_fus) +CD_wave*(pi*(self.d_f_outer/2)**2)/self.S
        
        return(CD0_fus, CDL_fus, CD_fus_sub, CD_fus_trans)
        
    def empennage_drag(self):
        
