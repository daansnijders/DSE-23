# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:23:05 2019

@author: Anique
"""

from math import * 
import numpy as np
import matplotlib.pyplot as plt

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
    def __init__(self,S,A,rho,rho_0,l_f,V_cruise,V_TO,mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes,CL_alpha,l_cockpit, l_cabin, l_tail,lambda_2_rad,lambda_4_rad,x_nlg,z_nlg,D_nlg,b_nlg,D_strutt_nlg,x_mlg,z_mlg,D_mlg,b_mlg,D_strutt_mlg,lambda_h_2_rad,lambda_v_2_rad, MAC_c, Cr_v, Ct_v, Cr_h, Ct_h, S_h, S_v, S_c, CL_alpha_h, de_da_h, i_h, alpha0L_h, A_h, CL_alpha_c, de_da_c, i_c, alpha0L_c, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c,S_ef, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap):
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
        self.lambda_4_rad   = lambda_4_rad
        self.x_nlg          = x_nlg
        self.z_nlg          = z_nlg
        self.b_nlg          = b_nlg
        self.D_nlg          = D_nlg
        self.D_strutt_nlg   = D_strutt_nlg
        self.x_mlg          = x_mlg
        self.z_mlg          = abs(z_mlg)
        self.b_mlg          = b_mlg
        self.D_mlg          = D_mlg
        self.D_strutt_mlg   = D_strutt_mlg
        self.lambda_h_2_rad = lambda_h_2_rad
        self.lambda_v_2_rad = lambda_v_2_rad
        self.MAC_c          = MAC_c
        self.Cr_h           = Cr_h                                 
        self.Ct_h           = Ct_h                                 
        self.Cr_v           = Cr_v                            
        self.Ct_v           = Ct_v
        self.S_h            = S_h
        self.S_v            = S_v
        self.S_c            = S_c
        self.CL_alpha_h     = CL_alpha_h
        self.de_da_h        = de_da_h
        self.i_h            = i_h
        self.alpha0L_h      = alpha0L_h
        self.A_h            = A_h
        self.CL_alpha_c     = CL_alpha_c
        self.de_da_c        = de_da_c
        self.i_c            = i_c
        self.alpha0L_c      = alpha0L_c
        self.A_c            = A_c
        self.l_fueltank     = l_fueltank
        self.d_fueltank     = d_fueltank
        self.l_nacel        = l_nacel
        self.d_nacel        = d_nacel
        self.i_n            = i_n
        self.delta_C_L_h    = delta_C_L_h
        self.delta_C_L_c    = delta_C_L_c
        self.S_ef           = S_ef
        self.SWF            = SWF
        self.SWF_LE         = SWF_LE
        self.Delta_C_L_flap = Delta_C_L_flap 
        self.b_flap         = b_flap
        self.b_slat         = b_slat

    def wing_drag(self):
        Re_f  = self.rho   * self.V_cruise * self.l_f / self.mu_37
        Re_f0 = self.rho_0 * self.V_TO     * self.l_f / self.mu_sl
        #These result in
        R_wf = 1.01         #(Figure 4.1) 
        
        cos_lambda_2_rad   = cos(self.lambda_2_rad)
        #This results in
        R_LS   = 1.21        #(Figure 4.2)
        R_LS_c = 1.21        #(Figure 4.2)
        
        L_prime = 1.2
        C_f_w = .00265
        t_c = 0.15
        
        """ Calculate chord at fuselage-wing intersection """
        x = self.taper_ratio * self.b/2 / (1 - self.taper_ratio)
        h2 = self.b / 2 + x
        h1 = h2 - self.d_f_outer/2
        c_fuselage_wing = h1/h2*self.Cr
        
        S_wet = 2*(2*((self.Ct + c_fuselage_wing) / 2 * (self.b/2 - self.d_f_outer/2)))
        
        C_D_0_w = R_wf * R_LS * C_f_w * (1 + L_prime * (t_c) + 100 * (t_c)**4) * S_wet/self.S
        
        """ C_D_L_w """
        r_LE = 0.687                  #Leading Edge radius
        RE_LER = self.rho * self.V_cruise * r_LE / self.mu_37
        R_par = RE_LER * 1/(tan(self.lambda_le_rad)) * sqrt(1 - (self.M*cos(self.lambda_le_rad))**2)
        R_par2 = self.A * self.taper_ratio / cos(self.lambda_le_rad) 
        #This results in 
        R = 0.95    #Figure 4.7
        
        beta = sqrt(1-self.M**2)
        c_l_alpha = np.rad2deg((1.3 + 0.5)/(7+9))
        k = c_l_alpha/(2*pi / beta)
        C_L_a_w = (2*pi*self.A)/(2 + ((self.A * beta / k)**2 * (1 + (tan(self.lambda_2_rad) / beta) ) + 4)**(1/2) )
                
        e = 1.1*(C_L_a_w / self.A)*(R * (C_L_a_w / self.A) + (1-R)*pi)
        
        C_L_w = 1.05 * self.CLdes
        C_D_L_w = C_L_w**2 / (pi * self.A * e)
        
        C_D_w_sub = C_D_0_w + C_D_L_w
        
        """ Transsonic drag """
        C_D_w_wave = 0.002      #Figure 4.11
        C_D_0_w_trans = C_D_0_w + C_D_w_wave
        
        CDL_CL2 = 0.025         #Figure 4.13
        C_D_L_w_trans = CDL_CL2 * C_L_w**2.
        
        C_D_w_trans = C_D_0_w_trans + C_D_L_w_trans
        
        return (C_D_w_sub, C_D_w_trans, C_D_0_w)
        
                        
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
        
        if self.l_cabin <= 23:
            CDL_fus = eta1*cdc*alpha**3*(Splf/self.S)
        else:
            CDL_fus = eta2*cdc*alpha**3*(Splf/self.S)
        
        CD_fus_sub = CD0_fus + CDL_fus
        
        
        CDf_fus = Cf_fus*(Swet_fus/self.S)
        CDp_fus = Cf_fus*(60/(self.l_f/self.d_f_outer)**3 + 0.0025*(self.l_f/self.d_f_outer))*Swet_fus/self.S
        CD_wave = 0.005     #Figure 4.22
        
        CD_fus_trans = Rwf*(CDf_fus + CDp_fus) +CD_wave*(pi*(self.d_f_outer/2)**2)/self.S
        
        return(CD_fus_sub, CD_fus_trans)

    def empennage_drag(self):
        #Subsonic for horizontal tail(h), vertical tail(v) and canard(c)
        #Zero lift drag calculations
        cos_lambda_v_2_rad   = cos(self.lambda_v_2_rad)
        cos_lambda_h_2_rad   = cos(self.lambda_h_2_rad)
        #This results in
        R_LS_h   = 1.21        #Figure 4.2
        R_LS_v   = 1.22        #Figure 4.2         
        R_LS_c   = 1.21        #Figure 4.2 
        
        C_h = (self.Cr_h + self.Ct_h)/2
        C_v = (self.Cr_v + self.Ct_v)/2
        Re_h_sub = self.rho_0   * self.V_TO * C_h / self.mu_sl
        Re_h_trans = self.rho   * self.V_cruise * C_h / self.mu_37
        Re_v_sub = self.rho_0 * self.V_TO     * C_v / self.mu_sl
        Re_v_trans = self.rho * self.V_cruise     * C_v / self.mu_37        
        Re_c_sub = self.rho_0 * self.V_TO     * self.MAC_c / self.mu_sl
        Re_c_trans = self.rho * self.V_cruise     * self.MAC_c / self.mu_37
        #This results in
        Cf_emp_h_sub = 0.0028       #Figure 4.3, config 1,2,3
        Cf_emp_h_trans = 0.0027     #Figure 4.3, config 1,2,3 
        Cf_emp_v_sub = 0.0026       #Figure 4.3, config 1,2,3 
        Cf_emp_v_trans = 0.0026     #Figure 4.3, config 1,2,3
        Cf_emp_c_sub_2 = 0.00325    #Figure 4.3, config 2
        Cf_emp_c_sub_3 = 0.0029     #Figure 4.3, config 3  
        Cf_emp_c_trans_2 = 0.00315  #Figure 4.3, config 2
        Cf_emp_c_trans_3 = 0.00285  #Figure 4.3, config 3 
        
        L_prime = 1.2           #Figure 4.4
        
        t_over_c = 0.15         #TO BE UPDATED AFTER AIRFOIL SELECTION EMPENNAGE AND CANARD
        
        Swet_h = 2*self.S_h
        Swet_v = 2*self.S_v
        Swet_c = 2*self.S_c
        
        
        CD0_h_tail_sub = R_LS_h*Cf_emp_h_sub*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_h/self.S
        CD0_v_tail_sub = R_LS_v*Cf_emp_v_sub*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_v/self.S
        
        if Re_c_sub == 0:
            CD0_c_sub = 0
        elif self.MAC_c <= 1.15:
            CD0_c_sub = R_LS_c*Cf_emp_c_sub_2*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_c/self.S
        else:
            CD0_c_sub = R_LS_c*Cf_emp_c_sub_3*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_c/self.S
        
        
        #Lift induced drag calculations
        e_h = 0.5
        e_c = 0.5
        CL0 = self.CL_alpha*(pi/180)*5.5
        alpha = (self.CLdes -CL0)/self.CL_alpha
        CDL_h_sub = ((self.CL_alpha_h*(alpha*(1-self.de_da_h)+self.i_h - self.alpha0L_h))**2)/(pi*self.A_h*e_h)*(self.S_h/self.S)
        
        if Re_c_sub == 0:
            CDL_c_sub = 0 
        else:
            CDL_c_sub = ((self.CL_alpha_c*(alpha*(1-self.de_da_c)+self.i_c - self.alpha0L_c))**2)/(pi*self.A_c*e_c)*(self.S_c/self.S)
        
        CDL_v_sub = 0
        
        CD_h_sub = CD0_h_tail_sub + CDL_h_sub
        CD_v_sub = CD0_v_tail_sub + CDL_v_sub
        CD_c_sub = CD0_c_sub + CDL_c_sub
        
        
        #Transonic for horizontal tail(h), vertical tail(v) and canard(c)
        #Zero lift drag calculations
        CD0_h_tail_trans = R_LS_h*Cf_emp_h_trans*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_h/self.S + 0.002*(self.S_h/self.S)
        CD0_v_tail_trans = R_LS_v*Cf_emp_v_trans*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_v/self.S + 0.002*(self.S_v/self.S)
        
        if Re_c_trans == 0:
            CD0_c_trans = 0
        elif self.MAC_c <= 1.15:
            CD0_c_trans = R_LS_c*Cf_emp_c_trans_2*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_c/self.S + 0.002*(self.S_c/self.S)
        else:
            CD0_c_trans = R_LS_c*Cf_emp_c_trans_3*(1+L_prime*t_over_c+100*(t_over_c)**4)*Swet_c/self.S + 0.002*(self.S_h/self.S)
        
        #Lift induced drag
        CDL_CL2 = 0.025     #Figure 4.13
        CL_h = self.CL_alpha_h*(alpha*(1-self.de_da_h)+self.i_h - self.alpha0L_h)
        CL_c = self.CL_alpha_c*(alpha*(1-self.de_da_c)+self.i_c - self.alpha0L_c)
        CDL_h_trans = CDL_CL2*CL_h**2
        CDL_c_trans = CDL_CL2*CL_c**2
        CDL_v_trans = 0
        
        CD_h_trans = CD0_h_tail_trans + CDL_h_trans
        CD_v_trans = CD0_v_tail_trans + CDL_v_trans
        CD_c_trans = CD0_c_trans + CDL_c_trans
        
        return(CD_h_sub, CD_v_sub, CD_c_sub, CD_h_trans, CD_v_trans, CD_c_trans)
        
    def nacelle_drag(self):
        Re_f  = self.rho   * self.V_cruise * self.l_nacel / self.mu_37
        Re_f0 = self.rho_0 * self.V_TO     * self.l_nacel / self.mu_sl
        #RE at service ceiling results in (same for all configurations)
        Rwf = 1.0               #Figure 4.1
        Cf_nacel = 0.00265    #Figure 4.3
        ratio = self.l_nacel/self.d_nacel
        Swet_nacel = pi*self.d_nacel*self.l_nacel*abs((1-2/ratio))**(2/3)*(1+1/(ratio)**2)
        
        CD0_nacel = 2*Rwf*Cf_nacel*(1+60/(self.l_nacel/self.d_nacel)**3 + 0.0025*(self.l_nacel/self.d_nacel))*Swet_nacel/self.S
        
        CL0 = self.CL_alpha*(pi/180)*5.5
        alpha = (self.CLdes - CL0)/self.CL_alpha
        alpha_n = alpha + self.i_n
        
        eta = 0.55          #Figure 4.19
        cdc = 1.2           #Figure 4.20
        Splf = self.d_nacel*self.l_nacel
        
        CDL_nacel = 2*eta*cdc*alpha_n**3*(Splf/self.S)
        
        CD_nacel_sub = CD0_nacel + CDL_nacel
        
        
        CDf_nacel = Cf_nacel*(Swet_nacel/self.S)
        CDp_nacel = Cf_nacel*(60/(self.l_nacel/self.d_nacel)**3 + 0.0025*(self.l_nacel/self.d_nacel))*Swet_nacel/self.S
        CD_wave = 0     #Figure 4.22
        
        CD_nacel_trans = 2*Rwf*(CDf_nacel + CDp_nacel) +CD_wave*(pi*(self.d_nacel/2)**2)/self.S
        
        return(CD_nacel_sub, CD_nacel_trans)
        
        
    def flaps_drag(self, C_D_0_W):
        """ Flaps """
        Cf_over_c = 0.35
        Delta_CDp_TO = 0.05
        Delta_CDp_land = 0.15
        Delta_CD_prof_TO = Delta_CDp_TO*np.cos(self.lambda_4_rad)*self.SWF/self.S
        Delta_CD_prof_land = Delta_CDp_land*np.cos(self.lambda_4_rad)*self.SWF/self.S
        
        Bfi_b = (self.d_f_outer + 1)/self.b
        Bfo_b = (self.d_f_outer + 1 + 2*self.b_flap)/self.b
        K = 0.22
        Delta_CD_i = K**2*self.Delta_C_L_flap**2*cos(self.lambda_4_rad)
        
        K_int = 0.4
        Delta_CD_int_TO = K_int*Delta_CD_prof_TO
        Delta_CD_int_land = K_int*Delta_CD_prof_land
        
        CD_flap_TO = Delta_CD_prof_TO + Delta_CD_i + Delta_CD_int_TO
        CD_flap_land = Delta_CD_prof_land + Delta_CD_i + Delta_CD_int_land
        
        """ Krueger """
        Cs_over_c = 1.1
        Delta_CDp_LE = C_D_0_W*Cs_over_c
        Delta_CD_prof_LE = Delta_CDp_LE*np.cos(self.lambda_4_rad)*self.SWF_LE/self.S
        
        Delta_CL_krug = 0.1
        Bfi_b_LE = (self.d_f_outer + 1)/self.b
        Bfo_b_LE = (self.d_f_outer + 1 + 2*self.b_slat)/self.b
        K_LE = 0.16
        Delta_CD_i_LE = K_LE**2*Delta_CL_krug**2*cos(self.lambda_4_rad)
        
        K_int_LE = 0.1
        Delta_CD_int_LE = K_int_LE*Delta_CD_prof_LE
        
        CD_slat = Delta_CD_prof_LE + Delta_CD_i_LE + Delta_CD_int_LE
        return(CD_flap_TO, CD_flap_land, CD_slat)
  
    def landinggear_drag(self):
        drag_par1 = self.x_nlg / self.D_nlg
        drag_par2 = self.z_nlg / self.D_nlg
        S_mlg = self.D_nlg * self.b_nlg * 2 + self.D_strutt_nlg * self.z_nlg
        #From this follows
        C_D_nlg = 0.5   #Figure 4.58
        
        S_mlg = self.D_mlg * self.b_mlg * 2 + self.D_strutt_mlg * self.z_mlg
        a = 2* self.b_mlg + self.D_strutt_mlg
        m = S_mlg / (a * self.z_mlg)
        #From this follows:
        C_D_mlg = 1.4   #Figure 4.59
        
        C_D_gear = (C_D_nlg ) * self.b_nlg * self.D_nlg  / self.S + 2*((C_D_mlg) * m / self.S)
        
        return (C_D_gear)
    
    
    def windshield_drag(self):
        delta_C_D_ws = 0.002        #Figure 4.68
        C_D_nosecone = 0.078        #Figure 4.68
        
        C_D_ws = delta_C_D_ws * (pi*(self.d_f_outer / 2)**2)/self.S
        
        return (C_D_ws)
    
    def store_drag(self):
        """ store_drag only applicable for configuration 3 """
        Re_f  = self.rho   * self.V_cruise * self.l_fueltank / self.mu_37
        Re_f0 = self.rho_0 * self.V_TO     * self.l_fueltank / self.mu_sl
        #RE at service ceiling results in (same for all configurations)
        Rwf = 0.95              #Figure 4.1
        Cf_fueltank = 0.003     #Figure 4.3
        ratio = self.l_fueltank/self.d_fueltank
        Swet_fueltank = pi*self.d_fueltank*self.l_fueltank*(1-2/ratio)**(2/3)*(1+1/(ratio)**2)
        
        CD0_fueltank = Rwf*Cf_fueltank*(1+60/(self.l_fueltank/self.d_fueltank)**3 + 0.0025*(self.l_fueltank/self.d_fueltank))*Swet_fueltank/self.S
        
        CL0 = self.CL_alpha*(pi/180)*5.5
        alpha = (self.CLdes -CL0)/self.CL_alpha
        
        eta = 0.66          #Figure 4.19
        cdc = 1.2           #Figure 4.20
        Splf = self.d_fueltank*self.l_fueltank
        
        CDL_fueltank = eta*cdc*alpha**3*(Splf/self.S)
        
        CD_fueltank_sub = 2* (CD0_fueltank + CDL_fueltank)
        
        
        CDf_fueltank = Cf_fueltank*(Swet_fueltank/self.S)
        CDp_fueltank = Cf_fueltank*(60/(self.l_fueltank/self.d_fueltank)**3 + 0.0025*(self.l_fueltank/self.d_fueltank))*Swet_fueltank/self.S
        CD_wave = 0.005     #Figure 4.22
        
        CD_fueltank_trans = 2*(Rwf*(CDf_fueltank + CDp_fueltank) +CD_wave*(pi*(self.d_fueltank/2)**2)/self.S)

        return(CD_fueltank_sub, CD_fueltank_trans)
        
        
    def trim_drag(self):
        
        if self.l_f <= 32 :
            delta_C_L_c = 0
        else: 
            if self.l_fueltank == 0:
                delta_C_L_c = 0
            else : 
                delta_C_L_c = self.delta_C_L_c
                
        e_h = 0.5           #if normal tail        
#        e_h = 0.75         #if T-tail
        e_c = 0.5
        
        if self.S_c == 0:
            delta_C_D_trim_lift = ((self.delta_C_L_h)**2 / (pi * self.A_h * e_h)) * self.S / self.S_h
        else: 
            delta_C_D_trim_lift = ((self.delta_C_L_h)**2 / (pi * self.A_h * e_h)) * self.S / self.S_h + ((delta_C_L_c)**2 / (pi * self.A_c * e_c)) * self.S / self.S_c
        
        #It follows that 
        delta_C_D_P_lambda_4_0 = 0.015       #Figure 4.44
        
        delta_C_D_trim_prof = delta_C_D_P_lambda_4_0 * cos(self.lambda_4_rad) * (self.S_ef / self.S_h)*(self.S_h / self.S)
        
        C_D_trim = delta_C_D_trim_lift + delta_C_D_trim_prof
        
        return (C_D_trim)
    
    def spoiler_drag(self):
        delta_C_D_sp = 0.63
        S_sp_i = [12,12]             #Area of the different spoilers
        C_D_sp = [delta_C_D_sp * (x / self.S) for x in S_sp_i]
        C_D_sp = sum(C_D_sp)
        
        return (C_D_sp)
    
    def surface_roughness_drag(self):
        k = 0.00000167
        ft_to_m=0.3048
        l_f_feet = self.l_f / ft_to_m        
        l_k = l_f_feet / k
        print (l_k)
        #From this, the cut-off Reynolds number follows
        Re_cutoff = 10**9 
        """ From this, it follows that there is no extra drag due to surface roughness """
        C_D_surface_roughness = 0
        return (C_D_surface_roughness)
    
    
class Lift:
<<<<<<< HEAD
    def __init__(self,S,A,rho,rho_0,l_f,V_cruise,M_cruise,V_TO,mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad, lambda_c_2_rad, i_c, S_h, S_c):
=======
    def __init__(self,S,A,rho,rho_0,l_f,V_cruise,M_cruise,V_TO,mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,i_w,wing_twist, A_h, A_c,lambda_h_2_rad, lambda_c_2_rad, i_c, S_h, S_c, i_h, x_le_MAC):
>>>>>>> 8d9e85f4b4f5213db9540796caf37facb1ca76d4
        self.S              = S
        self.A              = A
        self.rho            = rho
        self.rho_0          = rho_0
        self.l_f            = l_f
        self.V_cruise       = V_cruise
        self.M_cruise       = M_cruise
        self.V_TO           = V_TO
        self.mu_37          = mu_37
        self.MAC            = MAC
        self.Ct             = Ct
        self.Cr             = Cr
        self.b              = b
        self.taper_ratio    = taper_ratio
        self.d_f_outer      = d_f_outer
        self.lambda_le_rad  = lambda_le_rad
        self.lambda_4_rad   = lambda_4_rad
        self.lambda_2_rad   = lambda_2_rad
        self.alpha_0_l      = alpha_0_l
        self.C_l_alpha      = C_l_alpha
        self.alpha_C_l_max  = np.deg2rad(alpha_C_l_max)
        self.C_l_max        = C_l_max
<<<<<<< HEAD
        self.alpha_star_l   = alpha_star_l
        self.i_w            = i_w
        self.wing_twist     = wing_twist
=======
        self.i_w            = np.deg2rad(i_w)
        self.wing_twist     = np.deg2rad(wing_twist)
>>>>>>> 8d9e85f4b4f5213db9540796caf37facb1ca76d4
        self.A_h            = A_h
        self.A_c            = A_c
        self.lambda_h_2_rad = lambda_h_2_rad
        self.lambda_c_2_rad = lambda_c_2_rad
        self.i_c            = np.deg2rad(i_c)
        self.S_h            = S_h
        self.S_c            = S_c
        self.i_h            = i_h
        self.x_le_MAC       = x_le_MAC
        
        
    def Airfoil_lift_flaps(self):
        #Lift increase due to double slotted flaps
        #Wild guess for chord length flap one and two, together a little more than 0.35
        c1 = 0.2        #0.2*c
        c2 = 0.17       #0.17*c
        Phi_TE_upper = np.arctan(10*0.03)
        df1 = 40        #deg
        df2 = 10        #deg
        eta1 = 0.46     #Figure 8.20
        eta2 = 0.38     #Figure 8.20
        etat = 1.0      #Figure 8.22
        cldf1 = 0.0605   #Figure 8.21
        cldf2 = 0.054    #Figure 8.21
        
        #Wild guess for chord extension due to flaps
        c_a_prime = 1.13    #1.13*c
        c_prime = 1.20      #1.20*c

        delta_cl_flap = eta1*cldf1*df1*c_a_prime + eta2*etat*cldf2*df2*(1+(c_prime-c_a_prime))
        
        #Lift increase due to Krueger flaps
        cld = 0.0015    #Figure 8.26
        df = 60         #Wild guess
        c_prime_k = 1.1   
        
        delta_cl_krueger = cld*df*c_prime_k
        
         #lift curve slope
        c_prime_tot = 1.3 #c_prime and c_prime_k
        clalpha_flaps = c_prime_tot*self.C_l_alpha
        
        #Cl max increase due to TE and LE devices
        delclmax_base = 1.55    #Figure 8.31
        k1 = 1.2                #Figure 8.32
        k2 = 1.0                #Figure 8.33
        k3 = 1.0                #Figure 8.34
        delta_clmax_flap = delclmax_base*k1*k2*k3
        
        cldmax = 1.2            #Figure 8.35
        r_LE = 0.687            #Leading Edge radius
        eta_max = 0.82          #Figure 8.36
        df_rad = np.deg2rad(df) 
        
        delta_clmax_krueger = cldmax*eta_max*df_rad*c_prime_k
        
        return(delta_cl_flap, delta_cl_krueger, clalpha_flaps, delta_clmax_flap, delta_clmax_krueger)
       
    
    def Wing_lift(self):
        alpha = np.array([-1,0,1,3,5,7,9,11,13,15])  * pi/180 
        alpha_w = alpha - self.i_w
                
        alpha_0_l_m75_m30 = 0.033               #Figure 8.42
        
        delta_alpha_0 = (-0.343 - 0.398)/2      #Figure 8.41
        alpha_0_L_w = (self.alpha_0_l + delta_alpha_0 * self.wing_twist) * (alpha_0_l_m75_m30)
                
        C_l_alpha_M75 = self.C_l_alpha / sqrt(1 - self.M_cruise**2)
        beta = sqrt(1-self.M_cruise**2)      # Prandtl-Glauert compressibility correction factor
        k = C_l_alpha_M75 / (2*pi / beta)    # Constant dependent on the airfoil lift curve slope at M=0.75
        C_L_alpha_w = (2*pi*self.A)/(2 + sqrt(4 + (self.A*beta/k)**2 * (1+(tan(self.lambda_2_rad)**2)/beta**2)))
        print (C_L_alpha_w)
        
        """ Determining the spanwise lift distribution """
        L_b = np.array([(-0.293 + (-0.323 - -0.293)/2 * 1.5), (-0.204 + (-0.224 - -0.204)/2 * 1.5), (-0.012 + (-0.010 - -0.012)/2 * 1.5), (0.120 + (0.132 - 0.120)/2 * 1.5), (0.174 + (0.188 - 0.174)/2 * 1.5), (0.170 + (0.184 - 0.170)/2 * 1.5), (0.140 + (0.152 - 0.140)/2 * 1.5), (0.091 + (0.105 - 0.091)/2 * 1.5)])
        L_a = np.array([(1.392 + (1.409 - 1.392)/2 * 1.5), (1.294 + (1.299 - 1.294)/2 * 1.5), (1.150 + (1.148 - 1.150)/2 * 1.5), (0.956 + (0.947 - 0.956)/2 * 1.5), (0.710 + (0.704 - 0.710)/2 * 1.5), (0.536 + (0.541 - 0.536)/2 * 1.5), (0.403 + (0.410 - 0.403)/2 * 1.5), (0.283 + (0.295 - 0.283)/2 * 1.5)])
        y_b_2 = np.array([0,0.2,0.4,0.6,0.8,0.9,0.95,0.975])
        
        x = self.taper_ratio * self.b/2 / (1 - self.taper_ratio)
        h2 = self.b / 2 + x
        h1 = h2 - y_b_2 * self.b/2
        c = h1/h2*self.Cr
        
        eta = self.wing_twist
        a_0 = self.C_l_alpha
        E = 1.12
        a_e = a_0 / E
        
        C_L_w = C_L_alpha_w * (alpha_w - alpha_0_L_w)
        print (C_L_w)
        
        c_l_b1 = L_b * eta * a_e * self.S / c / self.b
        c_l_a1 = L_a * self.S / c / self.b
        
        c_l = []
        for i in range(len(alpha_w)):
            c_l_i = c_l_b1 + c_l_a1 * C_L_w[i]
            plt.plot(y_b_2, c_l_i)
            c_l.append(c_l_i)
            
        plt.show
        
        f = 0.997
        J = -0.38       #Figure 9  (Theory of Wing sections)
        e = 0.982       #Figure 10 (Theory of Wing sections)
        
        
        
        


#    def Wing_lift_flaps(self):
        
        
    
    def Airplane_lift(self, CL_alpha_w, alpha_0_L_w, CL_max_w, alpha_CL_max_w):
        beta = sqrt(1-self.M_cruise**2)
        CL_alpha_h = (2*pi*self.A_h)/(2 + sqrt(4+(self.A_h*beta/0.95)**2*(1+(np.tan(self.lambda_h_2_rad)**2)/beta**2)))
        CL_alpha_c = (2*pi*self.A_c)/(2 + sqrt(4+(self.A_c*beta/0.95)**2*(1+(np.tan(self.lambda_c_2_rad)**2)/beta**2)))
        etah = 0.9      #Source internet
        etac = 1.0
        
        Kwf = 1 + 0.025*(self.d_f_outer/self.b) - 0.25*(self.d_f_outer/self.b)**2
        CL_alpha_wf = Kwf*CL_alpha_w
        
        CL0 = (self.i_w - np.deg2rad(alpha_0_L_w))*CL_alpha_wf + CL_alpha_h*etah*(self.S_h/self.S)(self.i_h) + CL_alpha_c*etac*(self.S_c/self.S)(self.i_c)
        KA = (1/self.A) - 1/(1 + self.A**1.7)
        KL = (10 - 3*self.taper_ratio)/7
        l_h = 0.9*l_f - self.x_le_MAC - 0.25*self.MAC
        h_h = 0.75*self.d_f_outer
        Kh = (1-h_h/self.b)/(2*l_h/self.b)**(1/3)
        beta2 = sqrt(1-0**2)
        CL_alpha_w_M0 = (2*pi*self.A)/(2 + sqrt(4+(self.A*beta2/0.95)**2*(1+(np.tan(self.lambda_2_rad)**2)/beta2**2)))
        de_da = 4.44*((KA*KL*Kh*cos(self.lambda_4_rad)**0.5)**1.19)*(CL_alpha_w/CL_alpha_w_M0)
        de_da_c = 0.15  #Figure 8.67
        CL_alpha = CL_alpha_wf + CL_alpha_h*etah*(self.S_h/self.S)(1 - de_da) + CL_alpha_c*etac*(self.S_c/self.S)(1 + de_da_c)
        
        alpha_0_L = (-CL0)/CL_alpha
        
        delta_alpha_wc = np.deg2rad(3)          #rad
        alpha_CL_max = np.deg2rad(alpha_CL_max_w) - self.i_w - delta_alpha_wc
        
        CL_max = CL_max_w - CL_alpha_wf*delta_alpha_wc + CL_alpha_h*(self.S_h/self.S)*(alpha_CL_max*(1-de_da) + self.i_h) + CL_alpha_c*(self.S_c/self.S)*(alpha_CL_max*(1+de_da_c) + self.i_c)
        return(CL_alpha, alpha_0_L, CL_max)

#    def Airplane_lift_flaps(self, delta_CL_w): 
        