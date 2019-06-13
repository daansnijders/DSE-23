# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:50:12 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt
from inputs.constants import *
from inputs.concept_1 import *
from modules.Stability.cg_weight_config1 import x_cg_min1, x_cg_max1, x_le_MAC_range_perc, x_le_MAC_range
from modules.Stability.cg_weight_loadingdiagram import x_cg_min, x_cg_max,weight_pass
from modules.main_class2 import *

V_app = 70  #estimated by RB

class empennage:
    def __init__(self, config, x_ac, CL_a_h, CL_a_ah, de_da, S_h, l_h, S, c, Vh_V, x_le_MAC, Cm_ac, CL_ah, x_cg, CL_h):   
        self.config = config
        self.x_ac=x_ac #from nose in [m]
        self.CL_a_h = CL_a_h
        self.CL_a_ah = CL_a_ah
        self.de_da = de_da
        #self.S_h = S_h
        self.l_h_ac = l_h # can be removed
        self.S = S
        self.c = c
        self.Vh_V = Vh_V
        self.x_lemac = x_le_MAC
        self.Cm_ac = Cm_ac
        self.CL_ah = CL_ah
        self.x_cg = x_cg
        self.CL_h = CL_h

        #self.hortail_vol = self.S_h * self.l_h_ac / (self.S * self.c)
        
        self.plot_stability_tail(False)

    
    def calc_xnp(self):
        self.x_np = self.x_ac + (self.CL_a_h / self.CL_a_ah * (1-self.de_da) * self.hortail_vol * (self.Vh_V)**2)*self.c
        return self.x_np
    
    def calc_xcg(self):
        self.x_cg = self.calc_xnp() - 0.05*self.c
        return self.x_cg
    
    def calc_Cm(self):
        self.Cm_lesstail = self.Cm_ac + self.CL_ah * (self.x_cg-self.x_ac)/self.c
        self.Cm_tail = -self.CL_h * self.S_h * self.l_h_ac / (self.S * self.c) * (self.Vh_V)**2
        self.Cm = self.Cm_lesstail + self.Cm_tail
        return self.Cm    
    
    

    
    
    
    
    def plot_stability_tail(self, plot = True):
        """Stability excluding margin"""
        aa = 1/(self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h_ac*(self.Vh_V)**2)
        bb = (self.x_ac) / (self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h_ac*(self.Vh_V)**2)
        
        """Stability including margin"""
        dd = 1/(self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h_ac*(self.Vh_V)**2)
        ee = (self.x_ac - 0.05*self.c) / (self.CL_a_h/self.CL_a_ah*(1-self.de_da)*self.l_h_ac*(self.Vh_V)**2)
        
        """Controlability"""
        ff = self.CL_ah / (self.CL_h*self.l_h_ac*(self.Vh_V)**2)
        gg = (self.c*self.Cm_ac-self.CL_ah*self.x_ac)/(self.CL_h*self.l_h_ac*(self.Vh_V)**2)
        
        self.l = np.arange(x_le_MAC_range[0], (x_le_MAC_range[2]+self.c+0.01), 0.01)
        self.Sh_S1 = [] #stability xnp
        self.Sh_S2 = [] #stability xcg
        self.Sh_C1 = [] #controlability xac - Cmac/CL_ah
        
        for i in range (len(self.l)):
            self.Sh_S1.append(aa*self.l[i]-bb)
            self.Sh_S2.append(dd*self.l[i]-ee)
            self.Sh_C1.append(ff*self.l[i]+gg)
        
        
        
        def interpolate1(point1, point2):
            dydx = (point1[1] - point2[1]) / (point1[0] - point2[0])
            b = point1[1] - dydx * point1[0]
            return lambda y: (y-b) / dydx
        
        def interpolate2(point1, point2):
            dydx = (point1[1] - point2[1]) / (point1[0] - point2[0])
            b = point1[1] - dydx * point1[0]
            return lambda x: dydx * x + b
        
        f_min = interpolate1([x_cg_min1[0],x_le_MAC_range_perc[0]],[x_cg_min1[1],x_le_MAC_range_perc[1]])
        f_max = interpolate1([x_cg_max1[0],x_le_MAC_range_perc[0]],[x_cg_max1[1],x_le_MAC_range_perc[1]])
        
        f_S2 = interpolate2([self.l[0],self.Sh_S2[0]],[self.l[-1],self.Sh_S2[-1]])
        f_C1 = interpolate2([self.l[0],self.Sh_C1[0]],[self.l[-1],self.Sh_C1[-1]])
        
        diff = 1
        diff1 = 2
        y = 0.4215000000001215
        while abs(diff1) >= abs(diff) and abs(f_C1(f_min(y)) - f_S2(f_max(y))) > 0.000001:
            diff1 = f_S2(f_max(y))-f_C1(f_min(y))
            y += 0.0000001            
            if f_min(y) > x_cg_min1[1]:
                f_min = interpolate1([x_cg_min1[1],x_le_MAC_range_perc[1]],[x_cg_min1[2],x_le_MAC_range_perc[2]])
                f_max = interpolate1([x_cg_max1[1],x_le_MAC_range_perc[1]],[x_cg_max1[2],x_le_MAC_range_perc[2]])
            diff = f_S2(f_max(y))-f_C1(f_min(y))

        self.S_h_S = f_C1(f_min(y))
        x_le_MAC_l_f = y
        self.S_h = self.S_h_S * S
        self.x_le_MAC = x_le_MAC_l_f * l_f[0]

        
        if plot:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(x_cg_min1, x_le_MAC_range_perc)
            ax1.plot(x_cg_max1, x_le_MAC_range_perc)
            ax1.scatter(x_cg_min1, x_le_MAC_range_perc)
            ax1.scatter(x_cg_max1, x_le_MAC_range_perc)
            ax1.set(xlabel =  'x_cg', ylabel = 'x_le_MAC/l_f')
            
            ax1.scatter([f_min(y),f_max(y)],[y,y], color = 'b')
            ax1.plot([f_min(y),f_max(y)],[y,y], color = 'b')
    
            ax2 = ax1.twinx()
            ax2.plot(self.l, self.Sh_S1)
            ax2.plot(self.l, self.Sh_S2)
            ax2.plot(self.l, self.Sh_C1)
            ax2.set( ylim = [0,0.2565], ylabel = 'S_h/S')
            
            ax2.scatter([f_min(y),f_max(y)],[f_C1(f_min(y)),f_S2(f_max(y))], color = 'r')
            ax2.plot([f_min(y),f_max(y)],[f_C1(f_min(y)),f_S2(f_max(y))], color = 'r')

            
           
        
        # =============================================================================
        # Horizontal tail - NACA 63 010
        # =============================================================================

        
        self.A_h = 4.95                                                         # [-] aspect ratio horizontal tail
        self.taper_ratio_h = 0.39                                               # [-] taper ratio horizontal tail
        self.lambda_h_le_rad = np.deg2rad(34)                                   # [rad] leading edge sweep angle horizontal tail 
        self.t_c_h = 0.10                                                       # [-] tickness over chord ratio horizontal tail
        
        def get_x_h(l_f):
            return 0.9* l_f[0]
        
        def get_b(S, A):
            return np.sqrt(S*A)
        
        def get_Cr(S, taper_ratio, b):
            return 2*S/((1+taper_ratio)*b)
        
        def get_Ct(Cr, taper_ratio):
            return Cr * taper_ratio
        
        def get_lambda_4_rad_from_lambda_le(lambda_le_rad,Cr,b,taper_ratio):
            lambda_4_rad= np.arctan(np.tan(lambda_le_rad)+Cr/(2*b)*(taper_ratio-1))
            return lambda_4_rad
        
        def get_lambda_2_rad(lambda_4_rad,A,taper_ratio):
            return np.arctan(np.tan(lambda_4_rad)-1/A*(1-taper_ratio)/(1+taper_ratio))
        
        self.x_h = get_x_h(l_f)
        self.b_h = get_b(self.S_h, self.A_h)                                    # [m] span horizontal tail
        self.Cr_h = get_Cr(self.S_h, self.taper_ratio_h, self.b_h)              # [m] root chord length horizontal tail
        self.Ct_h = get_Ct(self.Cr_h, self.taper_ratio_h)                       # [m] tip chord length horizontal tail
        self.z_h = 0.75 * d_f_outer                                             # [m] height of the vertical tail
        self.l_h = 0.9*l_f[0] - self.x_le_MAC - 0.25*MAC                        # [m] distance 0.25mac-horizontal tail cg (still needs to be changed to class 2)
   
        self.lambda_h_4_rad = get_lambda_4_rad_from_lambda_le(self.lambda_h_le_rad,self.Cr_h,self.b_h,self.taper_ratio_h) # [rad] quarter chord sweep angle
        self.lambda_h_2_rad = get_lambda_2_rad(self.lambda_h_4_rad,self.A_h,self.taper_ratio_h) # [rad] half chord sweep angle

        # =============================================================================
        # Vertical tail - NACA 63 012
        # =============================================================================
                
        self.V_v = 0.1                                                          # [-] volume vertical tail
        self.A_v = 1.5                                                          # [-] aspect ratio vertical tail
        self.taper_ratio_v = 0.375                                              # [-] taper ratio vertical tail
        self.lambda_v_le_rad = np.deg2rad(40)                                   # [rad] leading edge sweep angle vertical tail
        self.t_c_v = 0.12                                                       # [-] tickness over chord ratio vertical tail
        self.x_v = self.x_h
        
        def get_S_v(S, b, x_cg, V_v, x_v):
            return [V_v*S* b / (x_v - x_cg[i]) for i in range(3)]
        

        self.S_v = min(get_S_v(S, b, x_cg, self.V_v, config1_cg.x_cg_vtail))    # [m^2] surface area vertical tail
        self.b_v = get_b(self.S_v, self.A_v)                                    # [m] span vertical tail
        self.Cr_v = get_Cr(self.S_v, self.taper_ratio_v, self.b_v)              # [m] root chord lengh vertical tail
        self.Ct_v = get_Ct(self.Cr_v, self.taper_ratio_v)                       # [m] tip chord length vertical tail
        
        self.lambda_v_4_rad = get_lambda_4_rad_from_lambda_le(self.lambda_v_le_rad,self.Cr_v,self.b_v,self.taper_ratio_v) # [rad] quarter chord sweep angle
        self.lambda_v_2_rad = get_lambda_2_rad(self.lambda_v_4_rad,self.A_v,self.taper_ratio_v) # [rad] half chord sweep angle
        
        
        
        # engine inoperative case
        N_e = thrust_max/2 * y_engine                                           # [N*m] moment caused by engine inoperative

        self.l_v = 0.9*l_f[0] - self.x_le_MAC - 0.25*MAC                        # [m] distance 0.25mac-vertical tail cg (still needs to be changed to class 2)

        C_y_max = 0.836                                                         # [-] maximum airfoil lift coefficient
        Y_v_max = C_y_max * 0.5*rho_0*V_app**2 * self.S_v                       # [N] force exerted by the vertical tail
        Y_v_req = N_e/self.l_v                                                  # [N] force required by the vertical tail
        C_y_req = Y_v_req/(0.5*rho_0*V_app**2*self.S_v)                         # [-] lift coefficient required vtail
        
        beta_max = 12.0                                                         # [deg] stall angle of the vertical tail
        beta_req = C_y_req / C_y_max * beta_max                                 # [deg] side-slip angle
        N_v_max = - Y_v_max * self.l_v                                          # [N*m] moment caused by the vertical tail

        assert N_e < -N_v_max                                                   # check if tail is capable enough
        
        return self.Sh_S1, self.Sh_S2, self.Sh_C1, self.S_h_S, self.x_le_MAC, self.S_h, self.z_h, self.l_v, self.l_h
    
    
    
e2 = empennage(1, (11.78+0.25*3.8), 3.82, 4.90, 0.3835, 21.72, 16., 93.5, 3.8, 1., 11.78, -0.3, 1.6, x_cg_max, -0.5838, )
    
    
    
    
    
    