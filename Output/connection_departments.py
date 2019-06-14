'inputs'
# from inputs.constants import *
# from inputs.performance_inputs import *
from inputs.concept_1 import *

'department modules'
from modules.Aerodynamics import *
from modules.sustainability.greenhousegasemissions import *
from modules.performance.class2_performance import *
from modules.sustainability.noise_defs import *


"""
AERODYNAMICS
"""
# initial guess for some values or needed from other departments, must be updated after iteration
i_w = 0                 #Angle of incidence wing in DEG
i_h = 0                 #Angle of incidence horinzontal tail in DEG
i_c = 0                 #Angle of incidence canard in DEG
i_n = 0                 #Angle of incidence nacelle in DEG
wing_twist = -3         #Wing twist in DEG
S_elev = 9.987716166863338      
Delta_C_L_flap = 2.4 - 1.5          #Cl_land - Cl_clean from performanca_input

'HLD'
config_HLD = HLD_class(Cl_land,Cl_clean,S,A,lambda_4_rad,taper_ratio,CL_alpha,lambda_le_rad,Cr,d_f_outer)
SWF, b_flap, SWF_LE, b_slat = config_HLD.HLD()


'Lift'
# initial values for lift
alpha_0_l = -5.4
alpha_star_l = 10        # from -7 to 3 deg
C_l_alpha = np.rad2deg(2/18)
alpha_C_l_max = np.rad2deg(10.75)
C_l_max = 1.58
#C_l_alpha_M75 = C_l_alpha / sqrt(1 - M_cruise**2)

#Set up different configurations
config1_Lift = Lift(S,A,rho,rho_0,l_f[0],V_cruise,M_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad[0], lambda_c_2_rad, i_c, S_h[0], S_c[0], i_h, x_le_MAC[0], b_flap, SWF)
config2_Lift = Lift(S,A,rho,rho_0,l_f[1],V_cruise,M_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad[1], lambda_c_2_rad, i_c, S_h[1], S_c[1], i_h, x_le_MAC[1], b_flap, SWF)
config3_Lift = Lift(S,A,rho,rho_0,l_f[2],V_cruise,M_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad,alpha_0_l,C_l_alpha,alpha_C_l_max,C_l_max,alpha_star_l,i_w,wing_twist, A_h, A_c,lambda_h_2_rad[2], lambda_c_2_rad, i_c, S_h[2], S_c[2], i_h, x_le_MAC[2], b_flap, SWF)

#Lift outputs
delta_cl_flap1, delta_cl_krueger1, clalpha_flaps1, delta_clmax_flap1, delta_clmax_krueger1 = config1_Lift.Airfoil_lift_flaps()
C_L_w1, CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1 = config1_Lift.Wing_lift()
delta_CL_w1, delta_CL_alpha_w1, delta_CL_max_w1 = config1_Lift.Wing_lift_flaps(delta_cl_flap1,CL_alpha_w1,C_l_alpha,(delta_clmax_flap1 + delta_clmax_krueger1),b_slat)
CL_alpha_h1, CL_alpha_c1, CL_alpha1, alpha_0_L1, CL_max1, de_da1, de_da_c1, alpha_CL_max1 = config1_Lift.Airplane_lift(CL_alpha_w1, alpha_0_L_w1, CL_max_w1, alpha_CL_max_w1)
delta_CL1, delta_CL_alpha1, delta_CL_max1 = config1_Lift.Airplane_lift_flaps(delta_CL_w1, CL_alpha_h1, CL_alpha_c1, delta_CL_alpha_w1, de_da1, delta_CL_max_w1)

delta_cl_flap2, delta_cl_krueger2, clalpha_flaps2, delta_clmax_flap2, delta_clmax_krueger2 = config2_Lift.Airfoil_lift_flaps()
C_L_w2, CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2 = config2_Lift.Wing_lift()
delta_CL_w2, delta_CL_alpha_w2, delta_CL_max_w2 = config2_Lift.Wing_lift_flaps(delta_cl_flap2,CL_alpha_w2,C_l_alpha,(delta_clmax_flap2 + delta_clmax_krueger2),b_slat)
CL_alpha_h2, CL_alpha_c2, CL_alpha2, alpha_0_L2, CL_max2, de_da2, de_da_c2, alpha_CL_max2 = config2_Lift.Airplane_lift(CL_alpha_w2, alpha_0_L_w2, CL_max_w2, alpha_CL_max_w2)
delta_CL2, delta_CL_alpha2, delta_CL_max2 = config2_Lift.Airplane_lift_flaps(delta_CL_w2, CL_alpha_h2, CL_alpha_c2, delta_CL_alpha_w2, de_da2, delta_CL_max_w2)

delta_cl_flap3, delta_cl_krueger3, clalpha_flaps3, delta_clmax_flap3, delta_clmax_krueger3 = config3_Lift.Airfoil_lift_flaps()
C_L_w3, CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3  = config3_Lift.Wing_lift()
delta_CL_w3, delta_CL_alpha_w3, delta_CL_max_w3 = config3_Lift.Wing_lift_flaps(delta_cl_flap3,CL_alpha_w3,C_l_alpha,(delta_clmax_flap3 + delta_clmax_krueger3),b_slat)
CL_alpha_h3, CL_alpha_c3, CL_alpha3, alpha_0_L3, CL_max3, de_da3, de_da_c3, alpha_CL_max3 = config3_Lift.Airplane_lift(CL_alpha_w3, alpha_0_L_w3, CL_max_w3, alpha_CL_max_w3)
delta_CL3, delta_CL_alpha3, delta_CL_max3 = config3_Lift.Airplane_lift_flaps(delta_CL_w3, CL_alpha_h3, CL_alpha_c3, delta_CL_alpha_w3, de_da3, delta_CL_max_w3)

'Drag'
#Values that must come from other departments
D_strutt_nlg = 0.15
D_strutt_mlg = 0.2
l_fueltank = 1.5
d_fueltank = 0.3

#Hard to calculate, not important for first iteration, maybe from Daan and Stijn
delta_C_L_h = 0                 #When i_h = 0 
delta_C_L_c = 0                 #When i_c = 0

#Set up different configurations
config1_Drag = Drag(S,A,rho,rho_0,l_f[0],V_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[0],CL_alpha1,l_cockpit, l_cabin[0], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, w_nlg, D_strutt_nlg, x_mlg[0], z_mlg, D_mlg, w_mlg, D_strutt_mlg, lambda_h_2_rad[0], lambda_v_2_rad[0], MAC_c[0], Cr_v[0], Ct_v[0], Cr_h[0], Ct_h[0], S_h[0], S_v[0], S_c[0], CL_alpha_h1, de_da1, i_h, alpha_0_L1, A_h, CL_alpha_c1, de_da_c1, i_c, alpha_0_L1, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
config2_Drag = Drag(S,A,rho,rho_0,l_f[1],V_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[1],CL_alpha2,l_cockpit, l_cabin[1], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, w_nlg, D_strutt_nlg, x_mlg[1], z_mlg, D_mlg, w_mlg, D_strutt_mlg, lambda_h_2_rad[1], lambda_v_2_rad[1], MAC_c[1], Cr_v[1], Ct_v[1], Cr_h[1], Ct_h[1], S_h[1], S_v[1], S_c[1], CL_alpha_h2, de_da2, i_h, alpha_0_L2, A_h, CL_alpha_c2, de_da_c2, i_c, alpha_0_L2, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)
config3_Drag = Drag(S,A,rho,rho_0,l_f[2],V_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,CLdes[2],CL_alpha3,l_cockpit, l_cabin[2], l_tail, lambda_2_rad, lambda_4_rad,x_nlg, z_nlg, D_nlg, w_nlg, D_strutt_nlg, x_mlg[2], z_mlg, D_mlg, w_mlg, D_strutt_mlg, lambda_h_2_rad[2], lambda_v_2_rad[2], MAC_c[2], Cr_v[2], Ct_v[2], Cr_h[2], Ct_h[2], S_h[2], S_v[2], S_c[2], CL_alpha_h3, de_da3, i_h, alpha_0_L3, A_h, CL_alpha_c3, de_da_c3, i_c, alpha_0_L3, A_c, l_fueltank, d_fueltank, delta_C_L_h, delta_C_L_c, S_elev, l_nacel, d_nacel, i_n, SWF, SWF_LE, Delta_C_L_flap, b_slat, b_flap)

#Drag outputs
CD_w_sub1, CD_w_trans1, CD_0_w1 = config1_Drag.wing_drag()
CD_fus_sub1, CD_fus_trans1 = config1_Drag.fuse_drag()
CD_h_sub1, CD_v_sub1, CD_c_sub1, CD_h_trans1, CD_v_trans1, CD_c_trans1 = config1_Drag.empennage_drag()
CD_nacel_sub1, CD_nacel_trans1 = config1_Drag.nacelle_drag()
CD_flap_TO1, CD_flap_land1, CD_slat1 = config1_Drag.flaps_drag(CD_0_w1)
CD_gear1 = config1_Drag.landinggear_drag()
CD_ws1 = config1_Drag.windshield_drag()
CD_store_sub1, C_D_store_trans1 = config1_Drag.store_drag()
CD_trim1 = config1_Drag.trim_drag()
CD_spoiler1 = config1_Drag.spoiler_drag()

CD_cruise1 = CD_w_trans1 + CD_fus_trans1 + CD_h_trans1 + CD_v_trans1 + CD_c_trans1 + CD_nacel_trans1 + CD_ws1 + C_D_store_trans1 + CD_trim1
CD_TO1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_TO1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1
CD_land1 = CD_w_sub1 + CD_fus_sub1 + CD_h_sub1 + CD_v_sub1 + CD_c_sub1 + CD_flap_land1 + CD_slat1 + CD_gear1 + CD_ws1 + CD_store_sub1 + CD_trim1

CD_w_sub2, CD_w_trans2, CD_0_w2 = config2_Drag.wing_drag()
CD_fus_sub2, CD_fus_trans2 = config2_Drag.fuse_drag()
CD_h_sub2, CD_v_sub2, CD_c_sub2, CD_h_trans2, CD_v_trans2, CD_c_trans2 = config2_Drag.empennage_drag()
CD_nacel_sub2, CD_nacel_trans2 = config2_Drag.nacelle_drag()
CD_flap_TO2, CD_flap_land2, CD_slat2 = config2_Drag.flaps_drag(CD_0_w2)
CD_gear2 = config2_Drag.landinggear_drag()
CD_ws2 = config2_Drag.windshield_drag()
CD_store_sub2, C_D_store_trans2 = config2_Drag.store_drag()
CD_trim2 = config2_Drag.trim_drag()
CD_spoiler2 = config2_Drag.spoiler_drag()

CD_cruise2 = CD_w_trans2 + CD_fus_trans2 + CD_h_trans2 + CD_v_trans2 + CD_c_trans2 + CD_nacel_trans2 + CD_ws2 + C_D_store_trans2 + CD_trim2
CD_TO2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_TO2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2
CD_land2 = CD_w_sub2 + CD_fus_sub2 + CD_h_sub2 + CD_v_sub2 + CD_c_sub2 + CD_flap_land2 + CD_slat2 + CD_gear2 + CD_ws2 + CD_store_sub2 + CD_trim2

CD_w_sub3, CD_w_trans3, CD_0_w3 = config3_Drag.wing_drag()
CD_fus_sub3, CD_fus_trans3 = config3_Drag.fuse_drag()
CD_h_sub3, CD_v_sub3, CD_c_sub3, CD_h_trans3, CD_v_trans3, CD_c_trans3 = config3_Drag.empennage_drag()
CD_nacel_sub3, CD_nacel_trans3 = config3_Drag.nacelle_drag()
CD_flap_TO3, CD_flap_land3, CD_slat3 = config3_Drag.flaps_drag(CD_0_w3)
CD_gear3 = config3_Drag.landinggear_drag()
CD_ws3 = config3_Drag.windshield_drag()
CD_store_sub3, C_D_store_trans3 = config3_Drag.store_drag()
CD_trim3 = config3_Drag.trim_drag()
CD_spoiler3 = config3_Drag.spoiler_drag()

CD_cruise3 = CD_w_trans3 + CD_fus_trans3 + CD_h_trans3 + CD_v_trans3 + CD_c_trans3 + CD_nacel_trans3 + CD_ws3 + C_D_store_trans3 + CD_trim3
CD_TO3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_TO3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3
CD_land3 = CD_w_sub3 + CD_fus_sub3 + CD_h_sub3 + CD_v_sub3 + CD_c_sub3 + CD_flap_land3 + CD_slat3 + CD_gear3 + CD_ws3 + CD_store_sub3 + CD_trim3


'Moments'
#Initial values
x_ref = 0.5                     #Own decision based on nothing
cl_des_airfoil = 0.651          #CL at zero angle of attack from JAVAfoil (RE = 17*10^6, 65-615)

#Set up different configurations
config1_Moment = Moment(S,A,rho,rho_0,l_f[0],V_cruise,M_cruise,V_TO[0],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad, t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger1, x_ref, cl_des_airfoil, wing_twist, y_MAC, C_L_w1, delta_CL_w1, SWF_LE, b_slat)
config2_Moment = Moment(S,A,rho,rho_0,l_f[1],V_cruise,M_cruise,V_TO[1],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad, t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger2, x_ref, cl_des_airfoil, wing_twist, y_MAC, C_L_w2, delta_CL_w2, SWF_LE, b_slat)
config3_Moment = Moment(S,A,rho,rho_0,l_f[2],V_cruise,M_cruise,V_TO[2],mu_37,mu_sl,MAC,Cr,Ct,b,taper_ratio,d_f_outer,lambda_le_rad,lambda_4_rad,lambda_2_rad, t_c, C_l_alpha, alpha_0_l, alpha_star_l,delta_cl_flap1,delta_cl_krueger3, x_ref, cl_des_airfoil, wing_twist, y_MAC, C_L_w3, delta_CL_w3, SWF_LE, b_slat)

cm_des_airfoil1, dcm_dcl_airfoil1 = config1_Moment.Airfoil_moment()
delta_cm_flap1, delta_cm_krueger1 = config1_Moment.Airfoil_moment_flaps(cm_des_airfoil1)
Cm0_w_sub1, Cm0_w_trans1, dCm_dCl_w1 = config1_Moment.Wing_moment()
delta_Cm_w_flaps1, delta_Cm_w_krueger1 = config1_Moment.Wing_moment_flaps(Cm0_w_sub1)

cm_des_airfoil2, dcm_dcl_airfoil2 = config2_Moment.Airfoil_moment()
delta_cm_flap2, delta_cm_krueger2 = config2_Moment.Airfoil_moment_flaps(cm_des_airfoil2)
Cm0_w_sub2, Cm0_w_trans2, dCm_dCl_w2 = config2_Moment.Wing_moment()
delta_Cm_w_flaps2, delta_Cm_w_krueger2 = config2_Moment.Wing_moment_flaps(Cm0_w_sub2)

cm_des_airfoil3, dcm_dcl_airfoil3 = config3_Moment.Airfoil_moment()
delta_cm_flap3, delta_cm_krueger3 = config3_Moment.Airfoil_moment_flaps(cm_des_airfoil3)
Cm0_w_sub3, Cm0_w_trans3, dCm_dCl_w3 = config3_Moment.Wing_moment()
delta_Cm_w_flaps3, delta_Cm_w_krueger3 = config3_Moment.Wing_moment_flaps(Cm0_w_sub3)

"""
CONTROL AND STABILITY
"""


"""
STRUCTURES
"""

"""
PERFORMANCE & PROPULSION
"""
'simulation accuracy inputs'
max_airport_altitude = 2000.  # [m]
altitude_resolution = 5  # number of different altitudes considered
mass_resolution = 20  # resolution of plotting mass vs take-off field length

'general inputs'
engine_failure = False
show_plots = True

# temporary values, will be removed as soon as aerodynamics inputs are ready
C_L_to = 1.9
C_L_la = 2.3
C_D_0 = 0.018117539865047032
C_D_to = 0.19542078310015343
C_D_la = 0.28017202214599246
C_L_cruise = 0.8
C_D_cruise = 0.05
lift_over_drag = LoverD[0]
aspect_ratio = A
oswald_efficiency_number = e

'analysis'
config1_Performance = Performance(C_L_to, C_L_la, C_L_cruise, C_D_0, C_D_to, C_D_la, C_D_cruise, S, OEW[0], MTOW[0], g,
                                  screen_height_to, screen_height_la, thrust_max, friction_coefficient_to,
                                  friction_coefficient_la, reverse_thrust_factor, engine_failure,
                                  thrust_setting_climb_out, thrust_setting_transition, M_payload[0], M_fuel[0],
                                  max_airport_altitude, altitude_resolution, mass_resolution, thrust_setting_climb,
                                  H_m, V_cruise, R[0], lift_over_drag, aspect_ratio,
                                  oswald_efficiency_number, correction_factor_to, show_plots)
#
# config2_Performance = Performance(C_L_to, C_L_la, C_L_cruise, C_D_0, C_D_to, C_D_la, C_D_cruise, S, OEW, MTOW, g,
#                                   screen_height_to, screen_height_la, thrust_max, friction_coefficient_to,
#                                   friction_coefficient_la, reverse_thrust_factor, engine_failure,
#                                   thrust_setting_climb_out, thrust_setting_transition, M_payload[1], M_fuel[1],
#                                   max_airport_altitude, altitude_resolution, mass_resolution, thrust_setting_climb,
#                                   H_m, V_cruise, R[1], lift_over_drag, aspect_ratio,
#                                   oswald_efficiency_number, correction_factor_to, show_plots)
#
# config3_Performance = Performance(C_L_to, C_L_la, C_L_cruise, C_D_0, C_D_to, C_D_la, C_D_cruise, S, OEW, MTOW, g,
#                                   screen_height_to, screen_height_la, thrust_max, friction_coefficient_to,
#                                   friction_coefficient_la, reverse_thrust_factor, engine_failure,
#                                   thrust_setting_climb_out, thrust_setting_transition, M_payload[2], M_fuel[2],
#                                   max_airport_altitude, altitude_resolution, mass_resolution, thrust_setting_climb,
#                                   H_m, V_cruise, R[2], lift_over_drag, aspect_ratio,
#                                   oswald_efficiency_number, correction_factor_to, show_plots)


"""
SUSTAINABILITY
"""
'Greenhouse Gas Emissions'
#CO_2_emissions=[get_CO2_emissions(M_fuel_burnt[i]) for i in range(3)]
#
#EI_NOx_phase=get_EI_NOX_fuelflow(fuel_flow)
#
#
#NOx_phase = [get_NOx_emissions_total(fuel_flow_phase,time_in_flight_phase,EI_NOx_phase)]
#
#
#Dp_Foo_NOx=[get_Dp_Foo_NOx_specific(NOx_total,T_phase)   ]
#
#NO_x_req=[get_NOx_reduction_CAEP(Dp_Foo_NOx_caep,Dp_Foo_NOx) ] #should come to 45%
#
#Fuel_per_pax=[get_CO2_per_passenger_per_km(M_fuel_burnt[i],N_pax[i],R[i]) for i in range(3)]


'Noise'
phi_observer=radians(1)
flap_deflection=np.radians(40)


b_flap=b/2*0.4
c_flap=0.35*Cr
S_flap=b_flap*c_flap

V_approach=64

r1,r2,r3,theta_1,theta_2,theta_3=simulate_flight_path(V_approach)

OSPL_dBA_tot_straight=EPNdB_calculations(r2,theta_2,phi_observer,V_approach, S_flap, b_flap,flap_deflection )
OSPL_dBA_tot_up=EPNdB_calculations(r1,theta_1,phi_observer,V_approach, S_flap, b_flap,flap_deflection )
OSPL_dBA_tot_down=EPNdB_calculations(r3,theta_3,phi_observer,V_approach, S_flap, b_flap,flap_deflection )