# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:56:29 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt

from inputs.concept_1 import *
from inputs.constants import *
from modules.Stability.cg_weight_config1 import *
from modules.Stability.cg_weight_loadingdiagram import *
from modules.Stability.control_surf_func import *
from modules.Stability.check_ground import *
from modules.Stability.empennage import *
from modules.testfile_aero import *
from modules.main_class2 import *



"""NEED FROM OTHER FILES"""
V_critical = 70                                                                 # [m/s] V1 speed/V_app
etah = 0.9                                                                      # [-] eta_h of aerodynamics
x_ac      = (x_le_MAC[0]+0.25*MAC)                                              # [m] x-location of the main wing ac
CL_a_h    = CL_alpha_h1                                                         # [-] CL_alpha_h
CL_a_ah   = CL_alpha_w1                                                         # [-] CL_alpha_(A-h)
de_da     = de_da                                                               # [-] downwash
Vh_V      = 1.                                                                  # [-] V_h/V velocity factors
Cm_ac     = -0.3                                        #TBD                    # [-] moment coefficient of main wing ac
CL_ah     = CL_max_w1                                                           # [-] CL_(A-h)
x_cg      = x_cg_max_flight1                                                    # [m] x-location of the most aft cg location for configuration 1 during flight
CL_h      = -0.8                                                                # [-] lift coefficient htail
CL_c      = 1.3                                         #zelf                   # [-] lift coefficient canard
CL_a_c    = CL_alpha_c2                                                         # [-] CL_alpha_canard
a_0       = alpha_0_l                                                           # [rad] zero lift angle of attack
i_h       = 0                                                                   # [rad] incidence angle htail
i_c       = 0                                                                   # [rad] incidence angle canard
CN_h_a    = CL_a_h                                                              # [-] C_N_h_alpha htail
CN_w_a    = CL_alpha_w1                                                         # [-] C_N_w_alpha main wing
CN_c_a    = CL_a_c                                                              # [-] C_N_c_alpha canard
CN_h_def  = 0.5                                         #zelf                   # [-] C_N_h_de elevator deflection
Vc_V      = 1.                                           #zelf                   # [-] V_c/V velocity factors
"""====================="""


# initialize class:
empennage1 = empennage(2, x_ac, CL_a_h, CL_a_ah, de_da, l_h[0], S, c, Vh_V, x_le_MAC[0], Cm_ac, CL_ah, x_cg, CL_h, CL_c, CL_a_c, a_0, i_h, i_c, CN_h_a, CN_w_a, CN_c_a, CN_h_def, Vc_V, V_critical)
empennage2 = empennage(3, x_ac, CL_a_h, CL_a_ah, de_da, l_h[0], S, c, Vh_V, x_le_MAC[0], Cm_ac, CL_ah, x_cg, CL_h, CL_c, CL_a_c, a_0, i_h, i_c, CN_h_a, CN_w_a, CN_c_a, CN_h_def, Vc_V, V_critical)


# outputs:
x_le_MAC        = empennage1.x_le_MAC_out                                       # [m] x-location of MAC main wing
x_le_MAC_l_f    = empennage1.x_le_MAC_l_f                                       # [-] xlemac over fuselage length
x_le_w          = get_le_wing(y_MAC,x_le_MAC, lambda_2_rad, MAC, Cr)            # [m] x-location of le main wing with updated lemac

S_h             = empennage1.S_h                                                # [m^2] surface area of htail
A_h             = empennage1.A_h                                                # [-] aspect ratio htail
b_h             = empennage1.b_h                                                # [m] span htail
Cr_h            = empennage1.Cr_h                                               # [m] root chord htail
Ct_h            = empennage1.Ct_h                                               # [m] tip chord htail
taper_ratio_h   = empennage1.taper_ratio_h                                      # [-] taper ratio htail
lambda_h_le_rad = empennage1.lambda_h_le_rad                                    # [rad] leading edge sweep htail
lambda_h_2_rad  = empennage1.lambda_h_2_rad                                     # [rad] half chord sweep htail
lambda_h_4_rad  = empennage1.lambda_h_4_rad                                     # [rad] quarter chord sweep htail
x_h             = empennage1.x_h                                                # [m] x-location of ac of the htail?
l_h             = empennage1.l_h                                                # [m] distance between c/4 on MAC of the main wing and horizontal tail

S_v             = empennage1.S_v                                                # [m^2] surface area of vtail
A_v             = empennage1.A_v                                                # [-] aspect ratio vtail
b_v             = empennage1.b_v                                                # [m] span vtail
Cr_v            = empennage1.Cr_v                                               # [m] root chord vtail
Ct_v            = empennage1.Ct_v                                               # [m] tip chord vtail
taper_ratio_v   = empennage1.taper_ratio_v                                      # [-] taper ratio vtail
lambda_v_le_rad = empennage1.lambda_v_le_rad                                    # [rad] leading edge sweep vtail
lambda_v_2_rad  = empennage1.lambda_v_2_rad                                     # [rad] half chord sweep vtail
lambda_v_4_rad  = empennage1.lambda_v_4_rad                                     # [rad] quarter chord sweep vtail
x_v             = empennage1.x_v                                                # [m] x-location of ac of the vtail?

taper_ratio_c2  = empennage1.taper_ratio_c                                      # [-] taper ratio canard
lambda_h_le_rad2= empennage1.lambda_h_le_rad                                    # [rad] leading edge sweep angle canard
t_c_c2          = empennage1.t_c_c                                              # [-] tickness over chord ratio canard   
Sc_S2           = empennage1.Sc_S                                               # [-] Ratio area canard (assumed for now)
S_c2            = empennage1.S_c                                                # [m^2] Surface area of the canard
A_c2            = empennage1.A_c                                                # [-] Aspect ratio of the canard
b_c2            = empennage1.b_c                                                # [m] span canard
Cr_c2           = empennage1.Cr_c                                               # [m] root chord length canard
Ct_c2           = empennage1.Ct_c                                               # [m] tip chord length canard
z_c2            = empennage1.z_c                                                # [m] veritcal height of the canard
l_c2            = empennage1.l_c                                                # [m] distance 0.25mac-wing to 0.25MAC canard    

taper_ratio_c3  = empennage2.taper_ratio_c                                      # [-] taper ratio canard
lambda_h_le_rad3= empennage2.lambda_h_le_rad                                    # [rad] leading edge sweep angle canard
t_c_c3          = empennage2.t_c_c                                              # [-] tickness over chord ratio canard   
Sc_S3           = empennage2.Sc_S                                               # [-] Ratio area canard (assumed for now)
S_c3            = empennage2.S_c                                                # [m^2] Surface area of the canard
A_c3            = empennage2.A_c                                                # [-] Aspect ratio of the canard
b_c3            = empennage2.b_c                                                # [m] span canard
Cr_c3           = empennage2.Cr_c                                               # [m] root chord length canard
Ct_c3           = empennage2.Ct_c                                               # [m] tip chord length canard
z_c3            = empennage2.z_c                                                # [m] veritcal height of the canard
l_c3            = empennage2.l_c                                                # [m] distance 0.25mac-wing to 0.25MAC canard   



# control surfaces: (inputs still need to be worked on...)
c_elev = get_c_elev(Cr_h, Ct_h, b_h)                                            # [m] chord length elevator
S_elev = get_S_elev(S_h)                                                        # [m^2] surface area elevator
b_elev = get_b_elev(S_elev,c_elev)                                              # [m] span elevator

c_rud = get_c_rud(Cr_v, Ct_v, b_v)                                              # [m] chord length rudder
S_rud = get_S_rud(S_v)                                                          # [m^2] surface area rudder
b_rud = get_b_rud(S_rud,c_rud)                                                  # [m] span rudder

c_ail = get_c_ail(Cr,Ct,b)                                                      # [m] chord length aileron
S_ail = get_S_ail(S)                                                            # [m^2] surface area aileron
b_ail = get_b_ail(b)                                                            # [m] span aileron

c_splr = get_c_splr(Cr, Ct, b)                                                  # [m] chord length spoiler
b_splr = get_b_splr(b)                                                          # [m] span spoiler


# Update cg's DIFFERENT CONFIG'S
# landing gear placement
x_mlg[0] = update_x_mlg(config1_cg.calc_z_cg(),theta_rad,beta_rad, x_cg_max_flight1, stroke,l_f[0]) # [m] x-location of the mlg
x_mlg[1] = max([x_mlg[0] + l_cutout, update_x_mlg(config2_cg.calc_z_cg(),theta_rad,beta_rad, x_cg_max_flight2, stroke,l_f[1])])
x_mlg[2] = max([x_mlg[0] + l_cutout, update_x_mlg(config3_cg.calc_z_cg(),theta_rad,beta_rad, x_cg_max_flight3, stroke,l_f[2])])

z_mlg = update_z_mlg(x_mlg[0],beta_rad,config1_cg.calc_x_cg(), config1_cg.calc_z_cg())                                    # [m] z-location of the mlg
#z_mlg=max(z_mlg)
#
#l_m = get_l_mw(x_mlg,x_cg)                                                      # [m] mlg distance from c.g
#l_n = get_l_nw(l_m,P_mw,N_mw,P_nw,N_nw)                                         # [m] nlg distance from c.g
#
#x_nlg = get_x_nlg(x_cg,l_n)                                                     # [m] x-location of nlg
#x_nlg=min(x_nlg)
#l_n=[x_cg[i]-x_nlg for i in range(3)]


#config1_ground      = check_ground(cg1_pass[0], cg2_pass[0], weight_pass[0], cg1_fuel[0], cg2_fuel[0], weight_fuel[0], x_nlg, x_mlg[0])     
#config2_ground      = check_ground(cg1_pass[1], cg2_pass[1], weight_pass[1], cg1_fuel[1], cg2_fuel[1], weight_fuel[1], x_nlg, x_mlg[1])     
#config3_ground      = check_ground(cg1_pass[2], cg2_pass[2], weight_pass[2], cg1_fuel[2], cg2_fuel[2], weight_fuel[2], x_nlg, x_mlg[2])     
#
#
#frac = np.ones((3,2))
#frac[0,0], frac[0,1], frac1 = config1_ground.check_equilibrium()
#frac[1,0], frac[1,1], frac2 = config2_ground.check_equilibrium()
#frac[2,0], frac[2,1], frac2 = config3_ground.check_equilibrium()



