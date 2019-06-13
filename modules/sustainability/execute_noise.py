# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:41:25 2019

@author: Lisa
"""
from modules.sustainability.noise_calc import *



r_observer=120
theta_observer=radians(90)
phi_observer=radians(1)




pe_2_flap =[get_effective_pressure_flap(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,L_flap,K_flap,a_flap,G_flap,flap_deflection) for f in centrefreq]


pe_2_slat =[ get_effective_pressure_slat(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,L_slat,K_slat,a_slat,G_slat)for f in centrefreq]

pe_2_wing = [get_effective_pressure_wing(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,L_wing,K_wing,a_wing,G_wing)for f in centrefreq]

pe_2_mlg =  [get_effective_pressure_lg(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,K_mlg,D_mlg,a_lg,G_mlg)for f in centrefreq]

pe_2_nlg =  [get_effective_pressure_lg(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,K_nlg,D_nlg,a_lg,G_nlg)for f in centrefreq]

pe_2_strut_main= [get_effective_pressure_strut(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,K_strut,L_strut_mlg,a_lg,G_mlg)for f in centrefreq]
pe_2_strut_nose= [get_effective_pressure_strut(f,rho_0,a_sl,M_TO,r_observer,theta_observer,phi_observer,K_strut,L_strut_nlg,a_lg,G_nlg)for f in centrefreq]


OSPL_dBA_flap=get_overall_sound_level_general(pe_2_flap,freq_delta,centrefreq)
OSPL_dBA_slat=get_overall_sound_level_general(pe_2_slat,freq_delta,centrefreq)
OSPL_dBA_wing=get_overall_sound_level_general(pe_2_wing,freq_delta,centrefreq)
OSPL_dBA_mlg=get_overall_sound_level_general(pe_2_mlg,freq_delta,centrefreq)
OSPL_dBA_nlg=get_overall_sound_level_general(pe_2_nlg,freq_delta,centrefreq)

if phi==0:
    OSPL_dBA_mlg_strut=0
    OSPL_dBA_nlg_strut=0
else:
    OSPL_dBA_mlg_strut=get_overall_sound_level_general(pe_2_strut_main,freq_delta,centrefreq)
    OSPL_dBA_nlg_strut=get_overall_sound_level_general(pe_2_strut_nose,freq_delta,centrefreq)


flap_fig=plot_pbl_vs_freq(pe_2_flap,centrefreq,freq_delta)
slat_fig=plot_pbl_vs_freq(pe_2_slat,centrefreq,freq_delta)
wing_fig=plot_pbl_vs_freq(pe_2_wing,centrefreq,freq_delta)
mlg_fig=plot_pbl_vs_freq(pe_2_mlg,centrefreq,freq_delta)
nlg_fig=plot_pbl_vs_freq(pe_2_nlg,centrefreq,freq_delta)


pe_tot=[pe_2_flap [i]+pe_2_slat[i] + pe_2_wing[i] +pe_2_mlg[i]+pe_2_nlg[i] +pe_2_strut_main[i]+ pe_2_strut_nose[i]for i in range(len(centrefreq))] 

OSPL_dBA_tot=get_overall_sound_level_general(pe_tot,freq_delta,centrefreq)



print('flap',OSPL_dBA_flap, 'dBA' )
print('slat',OSPL_dBA_slat, 'dBA'  )
print('wing',OSPL_dBA_wing, 'dBA'  )
print('mlg',OSPL_dBA_mlg , 'dBA' )
print('nlg',OSPL_dBA_nlg , 'dBA' )
print('m strut',OSPL_dBA_mlg_strut, 'dBA'  )
print('n strut',OSPL_dBA_nlg_strut, 'dBA' )

print('overall dBA',OSPL_dBA_tot, 'dBA' )