# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:50:12 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt
from inputs.constants import *
from inputs.concept_1 import *
from modules.Stability.Stability_runfile_empennage import x_cg_min1, x_cg_max1, x_le_MAC_range_perc, x_le_MAC_range
from modules.Stability.Stability_runfile import x_cg_min, x_cg_max
from modules.EXECUTE_FILE import *

V_app = 70  #estimated by RB

class empennage:
    def __init__(self, config, x_ac, CL_a_h, CL_a_ah, de_da, S_h, l_h, S, c, V_h, V, x_le_MAC, Cm_ac, CL_ah, x_cg, CL_h):   
        self.config = config
        self.x_ac=x_ac #from nose in [m]
        self.CL_a_h = CL_a_h
        self.CL_a_ah = CL_a_ah
        self.de_da = de_da
        self.S_h = S_h
        self.l_h = l_h
        self.S = S
        self.c = c
        self.V_h = V_h
        self.V = V
        self.x_lemac = x_le_MAC
        self.Cm_ac = Cm_ac
        self.CL_ah = CL_ah
        self.x_cg = x_cg
        self.CL_h = CL_h

        self.hortail_vol = self.S_h * self.l_h / (self.S * self.c)
        
        self.plot_stability_horitail(True)

    
    def calc_xnp(self):
        self.x_np = self.x_ac + (self.CL_a_h / self.CL_a_ah * (1-self.de_da) * self.hortail_vol * (self.V_h / self.V)**2)*self.c
        return self.x_np
    
    def calc_xcg(self):
        self.x_cg = self.calc_xnp() - 0.05*self.c
        return self.x_cg
    
    def calc_Cm(self):
        self.Cm_lesstail = self.Cm_ac + self.CL_ah * (self.x_cg-self.x_ac)/self.c
        self.Cm_tail = -self.CL_h * self.S_h * self.l_h / (self.S * self.c) * (self.V_h/self.V)**2
        self.Cm = self.Cm_lesstail + self.Cm_tail
        return self.Cm    
    
    def plot_stability_horitail(self, plot = True):
        """Stability excluding margin"""
        aa = 1/(self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h*(self.V_h/self.V)**2)
        bb = (self.x_ac) / (self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h*(self.V_h/self.V)**2)
        
        """Stability including margin"""
        dd = 1/(self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h*(self.V_h/self.V)**2)
        ee = (self.x_ac - 0.05*self.c) / (self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h*(self.V_h/self.V)**2)
        
        """Controlability"""
        ff = self.CL_ah / (self.CL_h*self.l_h*(self.V_h/self.V)**2)
        gg = (self.c*self.Cm_ac-self.CL_ah*self.x_ac)/(self.CL_h*self.l_h*(self.V_h/self.V)**2)
        
        self.l = np.arange(x_le_MAC_range[0], (x_le_MAC_range[2]+self.c+0.01), 0.01)
        self.Sh_S1 = [] #stability xnp
        self.Sh_S2 = [] #stability xcg
        self.Sh_C1 = [] #controlability xac - Cmac/CL_ah
        
        for i in range (len(self.l)):
            self.Sh_S1.append(aa*self.l[i]-bb)
            self.Sh_S2.append(dd*self.l[i]-ee)
            self.Sh_C1.append(ff*self.l[i]+gg)
        if plot:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(x_cg_min1, x_le_MAC_range_perc)
            ax1.plot(x_cg_max1, x_le_MAC_range_perc)
            ax1.scatter(x_cg_min1, x_le_MAC_range_perc)
            ax1.scatter(x_cg_max1, x_le_MAC_range_perc)
       
            
            ax2 = ax1.twinx()
            ax2.plot(self.l, self.Sh_S1)
            ax2.plot(self.l, self.Sh_S2)
            ax2.plot(self.l, self.Sh_C1)
            ax2.set( ylim =  (0,0.26))
            plt.show()
        
        """Put this on when iterating"""
        Sh_S = float(input("Input the optimal Sh/S ratio: "))
        x_le_MAC = float(input("Input the optimal Xlemac/Lf: ")) * l_f[0] 
        S_h = Sh_S * S
        #print (S_h)
        # =============================================================================
        # Horizontal tail - NACA 63 010
        # =============================================================================

        
        A_h = 4.95                                                              # [-] aspect ratio horizontal tail
        taper_ratio_h = 0.39                                                    # [-] taper ratio horizontal tail
        lambda_h_le_rad = np.deg2rad(34)                                        # [rad] leading edge sweep angle horizontal tail 
        t_c_h = 0.10                                                            # [-] tickness over chord ratio horizontal tail
        
        def get_x_h(l_f):
            return 0.9* l_f[0]
        
        def get_b(S, A):
            return np.sqrt(S*A)
        
        def get_Cr_h(S_h, taper_ratio_h, b_h):
            return 2*S_h/((1+taper_ratio_h)*b_h) 
        
        def get_Ct_h(Cr_h, taper_ratio_h):
            return Cr_h * taper_ratio_h
        
        x_h = get_x_h(l_f)
        b_h = get_b(S_h, A_h)                                                   # [m] span horizontal tail
        Cr_h = get_Cr_h(S_h, taper_ratio_h, b_h)                                # [m] root chord length horizontal tail
        Ct_h = get_Ct_h(Cr_h, taper_ratio_h)                                    # [m] tip chord length horizontal tail
        
#        print (b_h)
#        print (Cr_h)
#        print (Ct_h)
#        print (x_h)
        
        

        # =============================================================================
        # Vertical tail - NACA 63 012
        # =============================================================================
                
        self.V_v = 0.1                                                               # [-] volume vertical tail
        self.A_v = 1.5                                                               # [-] aspect ratio vertical tail
        self.taper_ratio_v = 0.375                                                   # [-] taper ratio vertical tail
        self.lambda_v_le_rad = np.deg2rad(40)                                        # [rad] leading edge sweep angle vertical tail
        self.t_c_v = 0.12                                                            # [-] tickness over chord ratio vertical tail
        
        def get_S_v(S, b, x_cg, V_v, x_v):
            return [V_v*S* b / (x_v - x_cg[i]) for i in range(3)]
        
        def get_Cr(S, taper_ratio, b):
            return 2*S/((1+taper_ratio)*b)
        
        def get_Ct(Cr, taper_ratio):
            return Cr * taper_ratio
        
        def get_lambda_4_rad_from_lambda_le(lambda_le_rad,Cr,b,taper_ratio):
            lambda_4_rad= np.arctan(np.tan(lambda_le_rad)+Cr/(2*b)*(taper_ratio-1))
            return lambda_4_rad
        
        def get_lambda_2_rad(lambda_4_rad,A,taper_ratio):
            return np.arctan(np.tan(lambda_4_rad)-1/A*(1-taper_ratio)/(1+taper_ratio))
        
        #print(get_S_v(S, b, x_cg, V_v, x_le_v))
        self.S_v = min(get_S_v(S, b, x_cg, self.V_v, config1_cg.x_cg_vtail))    # [m^2] surface area vertical tail
        self.b_v = get_b(self.S_v, self.A_v)                                    # [m] span vertical tail
        self.Cr_v = get_Cr(self.S_v, self.taper_ratio_v, self.b_v)              # [m] root chord lengh vertical tail
        self.Ct_v = get_Ct(self.Cr_v, self.taper_ratio_v)                       # [m] tip chord length vertical tail
        
        self.lambda_v_4_rad = get_lambda_4_rad_from_lambda_le(self.lambda_v_le_rad,self.Cr_v,self.b_v,self.taper_ratio_v) # [rad] quarter chord sweep angle
        self.lambda_v_2_rad = get_lambda_2_rad(self.lambda_v_4_rad,self.A_v,self.taper_ratio_v) # [rad] half chord sweep angle
        
        
        N_e = thrust_max/2 * y_engine                                           # [N*m] moment caused by engine inoperative
        l_v = 0.9*l_f[0] - x_le_MAC - 0.25*MAC                                  # [m] distance 0.25mac-vertical tail cg (still needs to be changed to class 2)

        C_y_max = 0.836                                                         # [-] maximum airfoil lift coefficient
        Y_v_max = C_y_max * 0.5*rho_0*V_app**2 * self.S_v                            # [N] force exerted by the vertical tail
        Y_v_req = N_e/l_v                                                       # [N] force required by the vertical tail
        C_y_req = Y_v_req/(0.5*rho_0*V_app**2*self.S_v)                              # [-] lift coefficient required vtail
        
        beta_max = 12.0                                                         # [deg] stall angle of the vertical tail
        beta_req = C_y_req / C_y_max * beta_max                                 # [deg] side-slip angle
        N_v_max = - Y_v_max * l_v                                               # [N*m] moment caused by the vertical tail

        assert N_e < -N_v_max                                                   # check if tail is capable enough
        


        
        return self.Sh_S1, self.Sh_S2, self.Sh_C1, Sh_S, x_le_MAC, S_h
    
e2 = empennage(1, (11.78+0.25*3.8), 3.82, 4.90, 0.3835, 21.72, 16., 93.5, 3.8, 1., 1., 11.78, -0.3, 1.6, x_cg_max, -0.5838, )

#r = e2.plot_stability_horitail(True)   

    






















#alpha0 = np.deg2rad(-5)
#C_m_ac_w = -0.3
#C_L_w_a = 2*np.pi
#C_D_w = 0.05
#
#class empennage(object):
#    def __init__(self,config,C_L_w,C_D_w,C_M_ac_w , aoa_rad,x_cg):
#        self.N_engine = T_req[config] * np.sin(i_e_rad)
#        self.T_engine = T_req[config] * np.cos(i_e_rad)
##        self.C_M_ac_w = C_M_ac_w
##        self.C_N_w = C_L_w * np.cos(aoa_rad) + C_D_w * np.sin(aoa_rad)
##        self.C_T_w = C_D_w * np.cos(aoa_rad) - C_L_w * np.sin(aoa_rad)
#        self.config = config
#        self.x_w = x_le_MAC + 0.25*MAC # 0.25?
#        self.z_w = t_c * MAC /2
#        self.x_cg = x_cg
#        
#    def calc_C_N(self): 
#        self.C_N = MTOW[self.config]/(0.5*rho*V_cruise**2 * S)*np.cos(theta_rad)
#        return self.C_N
#    
#    def calc_x_h(self):
#        if self.config == 1:
#            A = (self.calc_C_N() - self.C_N_w - self.N_engine/(0.5*rho*V_cruise**2*S)*np.cos(i_e_rad) - self.T_engine/(0.5*rho*V_cruise**2*S)*np.sin(i_e_rad))
#            
#            x_h = self.x_cg + self.C_M_ac_w *MAC/A + self.C_N_w * (self.x_cg - self.x_w)/A - self.C_T_w * (z_cg - self.z_w) + (self.T_engine * np.cos(i_e_rad) - self.N_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (z_cg - z_engine) + (self.N_engine * np.cos(i_e_rad) + self.T_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (self.x_cg - x_engine)
#        return x_h
#    
#    def findeq(self):
#        aoa_rad = np.linspace(np.deg2rad(0),np.deg2rad(5),1000)
#        y1 = MTOW[self.config]*g / (0.5*rho*V_cruise**2*S) * np.cos(aoa_rad) - C_L_w_a * (aoa_rad - alpha0) * np.cos(aoa_rad) - C_D_w * np.sin(aoa_rad) - self.N_engine / (0.5*rho*V_cruise**2*S) * np.cos(i_e_rad) - self.T_engine/(0.5*rho*V_cruise**2*S)*np.sin(i_e_rad)
#        y2 = (-C_m_ac_w - (C_L_w_a * (aoa_rad - alpha0) * np.cos(aoa_rad) - C_D_w * np.sin(aoa_rad)) * (x_cg[self.config] - self.x_w[self.config])/MAC + (C_D_w * np.cos(aoa_rad) - C_L_w_a * (aoa_rad - alpha0) * np.sin(aoa_rad)) * (z_cg - self.z_w) / MAC - (self.T_engine * np.cos(i_e_rad) - self.N_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (z_cg - z_engine) + (self.N_engine * np.cos(i_e_rad) + self.T_engine * np.sin(i_e_rad)) / (0.5*rho*V_cruise**2*S*MAC) * (self.x_cg - x_engine[self.config])) * MAC / (self.x_cg - x_le_h[self.config])
#    
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#        ax.plot(np.rad2deg(aoa_rad),y1)
#        ax.plot(np.rad2deg(aoa_rad),y2)
#        
#        start = y1[0] < y2[0]
#        condition = start
#        i = 1
#        while condition == start:
#            condition = y1[i] < y2[i]
#            i += 1
#        
#        aoa_cruise = aoa_rad[i]
#        y_cruise = np.sum([y1[i],y2[i]])/2
#        
#        
#        return np.rad2deg(aoa_cruise),y_cruise


"""
Values needed
M = 0
M_ac_w = sybren
M_ac_h = 0
*M_ac_c
C_N_w = sybren
C_N_h = sybren
*N_c 
T_w = tbc wrt N_w
T_h = tbc wrt N_h
*T_c
x_cg = get_cg
z_cg = get_cg
*x_w
z_w = around 0 + 0.5 * t
x_h = 0.9 of l_f?
z_h = 0.8 * d_f_outer
*x_c
*z_c
N_e = np.sin(i_e) * thrust
T_e = np.cos(i_e) * thrust
i_e = 0
x_e = get_cg
z_e = get_cg
V_cruise = concept 1
V_h = ?
V_c = ?
theta/aoa = aerodynamic guys -> alpha cruise
we need cl/alpha in order to determine the lift produced.
"""