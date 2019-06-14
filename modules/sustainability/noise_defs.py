# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:53:00 2019

@author: Lisa
"""
from math import *
from inputs.concept_1 import pi, S, b, N_mw, N_nw, L_strut_mlg, L_strut_nlg
from inputs.constants import mu_sl, rho_0,a_sl,D_mlg,D_nlg
import matplotlib.pyplot as plt
#noise definitions


'set up the effective pressure of the different airframe components'
#f will be a list which will be an 1./third octave band or complete octave band 

def simulate_flight_path(V_approach):
    landing_angle=3*pi/180
    straight_distance=5*V_approach
    h_2=2300*tan(landing_angle)
    h_1=tan(landing_angle)*straight_distance+h_2
    h_3=-tan(landing_angle)*straight_distance+h_2
    
    r2=h_2
    r3=h_3/cos(landing_angle)
    r1=h_1/cos(landing_angle)
    theta_1=87*pi/180
    theta_2=pi/2
    theta_3=93*pi/180
    
    return r1,r2,r3,theta_1,theta_2,theta_3

def get_constants_wing(M):
    G_wing=0.37*S/b**2*((rho_0*M*a_sl*S)/(mu_sl*b))**-0.2
    L_wing=G_wing*b
    K_wing=4.464E-5 
    a_wing=5
    return G_wing,L_wing,K_wing,a_wing

def get_constants_slats(M):
    G_slat=0.37*S/b**2*((rho_0*M*a_sl*S)/(mu_sl*b))**-0.2
    L_slat=G_slat*b
    K_slat=4.464E-5 
    a_slat=5
    return G_slat,L_slat,K_slat,a_slat


def get_constants_flaps(M,flap_deflection,S_flap,b_flap):
    G_flap=S_flap/b**2*(sin(flap_deflection))**2
    L_flap=S_flap/b_flap
    K_flap=2.787*10**-4
    a_flap=6
    return G_flap,L_flap,K_flap,a_flap

def get_constants_landinggear(N_mw,N_nw,D_mlg,D_nlg):
    G_mlg=N_mw*(D_mlg/b)**2
    G_nlg=N_nw*(D_nlg/b)**2
    K_mlg=3.414E-4  #4 wheels
    K_nlg=4.349E-4 #1 or 2 wheels (nose)
    K_strut=2.753E-4 
    a_lg=6
    return G_mlg, G_nlg, K_mlg,K_nlg,K_strut,a_lg

def get_octave_frequency_bands():
    number_of_bands=14
    bandnumbers=list(range(1,number_of_bands+1))
    centrefreq=[10**(3*bandnumbers[i]/10) for i in range(len(bandnumbers))]
    lowerfreq=[2**(-0.5)*centrefreq[i] for i in range(len(bandnumbers))]
    upperfreq=[2**(0.5)*centrefreq[i] for i in range(len(bandnumbers))]
    freq_delta=[upperfreq[i]-lowerfreq[i] for i in range(len(bandnumbers))]
    return centrefreq,freq_delta

'strouhal numbers'
def get_strouhal_number(f,L,M,theta,a):
    S=f*L*(1-M*cos(theta))/(M*a)
    return S
def get_strouhal_number_lg_and_strut(f,d_w,M,a):
    S=f*d_w/(M*a) # 
    return S


'get acoustic powers'
def get_noise_acoustic_power(K,M,a_const,G,rho,a):
    P=K*M**a_const*G*(rho*a**3*b**2)
    return P
def get_noise_acoustic_power_landinggear(K,M,N_w,d_w):
    return K*M**6*N_w*d_w**2

def get_noise_acoustic_power_strut(K,M,d_w,L_strut):
    return K*M**6*d_w*L_strut


'spectrical functions'
def get_spectrical_function_wing(S):
    F=0.613*(10*S)**4*((10*S)**1.5+0.5)**-4
    return F
def get_spectrical_function_flap(S):
    if S<2:
        F=0.048*S
    if S>=2 and S<=20:
        F=0.1406*S**-0.55
    else:
        F=216.49*S**-3
    return F
def get_spectrical_function_slat(S):
    F=0.613*(10*S)**4*((10*S)**1.5+0.5)**-4 + 0.613*(2.19*S)**4*((2.19*S)**1.5+0.5)**-4
    return F
def get_spectrical_function_landinggear(S): #this is for two wheels
    F=2*(13.59*S**2*(S**2+12.5)**-2.25)
    return F






'get directivity functions'
def get_directivity_function_wing(theta,phi):  #same for slat
    D=4*(cos(phi))**2*(cos(theta/2))**2
    return D

def get_directivity_function_flap(flap_deflection,theta,phi):
    D=3*(sin(flap_deflection)*cos(theta)+cos(flap_deflection)*sin(theta)*cos(phi))**2
    return D

def get_directivity_function_landinggear(theta):
    D=3/2*(sin(theta))**2
    return D

def get_directivity_function_strut(theta,phi):
    D=3*(sin(theta))**2*(sin(phi))**2
    return D






'this is to get the effective pressure level (squared)'

def get_effective_pressure_flap(f,rho,a,M,r_observer,theta,phi,L,K,a_const,G,flap_deflection):
    Strouhal=get_strouhal_number(f,L,M,theta,a)
    F=get_spectrical_function_flap(Strouhal)
    P=get_noise_acoustic_power(K,M,a_const,G,rho,a)
    D=get_directivity_function_flap(flap_deflection,theta,phi)
    
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared

def get_effective_pressure_slat(f,rho,a,M,r_observer,theta,phi,L,K,a_const,G):
    Strouhal=get_strouhal_number(f,L,M,theta,a)
    F=get_spectrical_function_slat(Strouhal)
    P=get_noise_acoustic_power(K,M,a_const,G,rho,a)
    D=get_directivity_function_wing(theta,phi)
    
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared

def get_effective_pressure_wing(f,rho,a,M,r_observer,theta,phi,L,K,a_const,G):
    Strouhal=get_strouhal_number(f,L,M,theta,a)
    F=get_spectrical_function_wing(Strouhal)
    P=get_noise_acoustic_power(K,M,a_const,G,rho,a)
    D=get_directivity_function_wing(theta,phi)
    
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared

def get_effective_pressure_lg(f,rho,a,M,r_observer,theta,phi,K,d_w,a_const,G):
    Strouhal=get_strouhal_number(f,d_w,M,theta,a)
    F=get_spectrical_function_landinggear(Strouhal)
    P=get_noise_acoustic_power(K,M,a_const,G,rho,a)
    D=get_directivity_function_landinggear(theta)
    
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared

def get_effective_pressure_strut(f,rho,a,M,r_observer,theta,phi,K,L_strut,a_const,G):
    Strouhal=get_strouhal_number(f,L_strut ,M,theta,a)
    F=get_spectrical_function_landinggear(Strouhal)
    P=get_noise_acoustic_power(K,M,a_const,G,rho,a)
    D=get_directivity_function_strut(theta,phi)
    
    p_e_squared=rho*a*P*D*F/(4*pi*r_observer**2*(1-M*cos(theta))**4) #(pa^2)
    return p_e_squared



'get the sound pressure level'
def get_sound_pressure_level(p_e_squared):   #[dB]
    SPL=10*log10(p_e_squared/(2E-5)**2)
    return SPL

'correct for the third octave bands'
def correction_for_thirdbandwidth(SPL,freq_delta):
    SPL_bandwidth=SPL+10*log10(freq_delta)
    return SPL_bandwidth

def atmospheric_absorption(SPL_bandwidth,r,centrefreq):
    absorp_coeff=[0,0,0,0,0,0.122,0.440,1.31,2.73,4.67,9.9,29.7,106,364]
    SPL_bandwidth=[SPL_bandwidth[i]-absorp_coeff[i]*r/1000  for i in range(len(centrefreq))]
    return SPL_bandwidth



'perform A weightening '  
def A_weighting_correction(f):
    delta_L_A=-145.528+98.262*log10(f)-19.509*(log10(f))**2+0.975*(log10(f))**3
    return delta_L_A
 

def A_weighted_sound_level(f,delta_L_A,SPL):
    L_A_perfreq=SPL+delta_L_A 
    return L_A_perfreq

def A_weighted_overall_sound_level(centrefreq,L_A_perfreq):
    exponent=[10**(L_A_perfreq[i]/10) for i in range(len(centrefreq))]
    sum_exponent=sum(exponent)
    OSPL_A=10*log10(sum_exponent)

    return OSPL_A


'perform all in one'
def get_overall_sound_level_general(p_e_squared,freq_delta,centrefreq,r):
    SPL         =[get_sound_pressure_level(p_e_squared[i]) for i in range(len(centrefreq))]
    PBL         =[correction_for_thirdbandwidth(SPL[i],freq_delta[i]) for i in range(len(centrefreq))]
    PBL_absor         = atmospheric_absorption(PBL,r,centrefreq)
    
    delta_L_A   = [A_weighting_correction(centrefreq[i]) for i in range(len(centrefreq))]
    PBL_a = [A_weighted_sound_level(centrefreq,delta_L_A[i],PBL_absor[i])for i in range(len(centrefreq))]
    
    OSPL_A      =A_weighted_overall_sound_level(centrefreq,PBL_a)
    
    
    plt.figure()
    plt.plot(centrefreq,PBL_absor,'-o',label=' PBL absoption')
    plt.plot(centrefreq,PBL,'-o',label='PBL')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('1/1 octave PBL [dB(A)]')
    plt.xscale('log')
    plt.title('1/1 octave Pressure Band Level ')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(centrefreq,delta_L_A)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('$\\Delta$ L_A')
    plt.xscale('log')
    plt.xlim([100,10**4])
    plt.show()
    
    return OSPL_A



def plot_pbl_vs_freq(p_e_squared,centrefreq,freq_delta):
    SPL         =[get_sound_pressure_level(p_e_squared[i]) for i in range(len(centrefreq))]
    PBL         =[SPL[i]+10*log10(freq_delta[i]) for i in range(len(centrefreq))]
    delta_L_A   = [A_weighting_correction(centrefreq[i]) for i in range(len(centrefreq))]
    PBL_a = [A_weighted_sound_level(centrefreq,delta_L_A[i],PBL[i])for i in range(len(centrefreq))]
    
    
    

    plt.figure()
    plt.plot(centrefreq,PBL, '-o',label='Unweighted PBL')
    plt.plot(centrefreq,PBL_a,'-o',label='A-weighted PBL')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('1/3 octave PBL [dB(A)]')
    plt.xscale('log')
    plt.xlim(10,10**4)
    plt.title('1/3 octave Pressure Band Level A-weightening comparison')
    plt.legend()
    plt.show()




def get_the_limit_for_lateral(MTOW):
    EPNdB_limit=94+log10(35000)+(9)/(log10(400000)-log10(35000))*(log10(MTOW)-log10(35000))
    return EPNdB_limit


'redundant definitions for now '
def get_overall_sound_pressure_level(SPL,centrefreq):
    OSPL=10*log10(sum(10**(SPL/10)))
    return OSPL
    
def get_perceived_noise_level(N):
    L_pn=40+33.3*log10(N)
    return L_pn

def get_effective_perceived_noise_level(L_pn,delta_t):
    L_epn=10*log10(delta_t/10*sum(10**(L_pn/10)))
    return L_epn


def doppler_effect_moving_away(f,V_approach):
    f=f/(1+V_approach/a_sl)
    return f

def doppler_effect_moving_towards(f,V_approach):
    f=f/(1-V_approach/a_sl)
    return f 

def EPNdB_calculations(r_observer,theta_observer,phi_observer,V_approach, S_flap, b_flap,flap_deflection ):
    
    
    
    centrefreq,freq_delta=get_octave_frequency_bands()
    
    M=V_approach/a_sl
    
    if r_observer>121.:
        centrefreq=[doppler_effect_moving_towards(f,V_approach) for f in centrefreq]
    else:
        centrefreq=[doppler_effect_moving_away(f,V_approach) for f in centrefreq]
        
    

    G_wing,L_wing,K_wing,a_wing=get_constants_wing(M)
    G_slat,L_slat,K_slat,a_slat=get_constants_slats(M)
    G_flap,L_flap,K_flap,a_flap=get_constants_flaps(M,flap_deflection,S_flap,b_flap)
    G_mlg, G_nlg, K_mlg,K_nlg,K_strut,a_lg=get_constants_landinggear(N_mw,N_nw,D_mlg,D_nlg)


    pe_2_flap =[get_effective_pressure_flap(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,L_flap,K_flap,a_flap,G_flap,flap_deflection) for f in centrefreq]


    pe_2_slat =[ get_effective_pressure_slat(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,L_slat,K_slat,a_slat,G_slat)for f in centrefreq]

    pe_2_wing = [get_effective_pressure_wing(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,L_wing,K_wing,a_wing,G_wing)for f in centrefreq]

    pe_2_mlg =  [get_effective_pressure_lg(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,K_mlg,D_mlg,a_lg,G_mlg)for f in centrefreq]

    pe_2_nlg =  [get_effective_pressure_lg(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,K_nlg,D_nlg,a_lg,G_nlg)for f in centrefreq]

    pe_2_strut_main= [get_effective_pressure_strut(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,K_strut,L_strut_mlg,a_lg,G_mlg)for f in centrefreq]
    pe_2_strut_nose= [get_effective_pressure_strut(f,rho_0,a_sl,M,r_observer,theta_observer,phi_observer,K_strut,L_strut_nlg,a_lg,G_nlg)for f in centrefreq]


    OSPL_dBA_flap=get_overall_sound_level_general(pe_2_flap,freq_delta,centrefreq,r_observer)
    OSPL_dBA_slat=get_overall_sound_level_general(pe_2_slat,freq_delta,centrefreq,r_observer)
    OSPL_dBA_wing=get_overall_sound_level_general(pe_2_wing,freq_delta,centrefreq,r_observer)
    OSPL_dBA_mlg=get_overall_sound_level_general(pe_2_mlg,freq_delta,centrefreq,r_observer)
    OSPL_dBA_nlg=get_overall_sound_level_general(pe_2_nlg,freq_delta,centrefreq,r_observer)

    if phi_observer==0:
        OSPL_dBA_mlg_strut=0
        OSPL_dBA_nlg_strut=0
    else:
        OSPL_dBA_mlg_strut=get_overall_sound_level_general(pe_2_strut_main,freq_delta,centrefreq,r_observer)
        OSPL_dBA_nlg_strut=get_overall_sound_level_general(pe_2_strut_nose,freq_delta,centrefreq,r_observer)


    pe_tot=[pe_2_flap [i]+pe_2_slat[i] + pe_2_wing[i] +pe_2_mlg[i]+pe_2_nlg[i] +pe_2_strut_main[i]+ pe_2_strut_nose[i]for i in range(len(centrefreq))] 

    OSPL_dBA_tot=get_overall_sound_level_general(pe_tot,freq_delta,centrefreq,r_observer)



    print('flap',OSPL_dBA_flap, 'dBA' )
    print('slat',OSPL_dBA_slat, 'dBA'  )
    print('wing',OSPL_dBA_wing, 'dBA'  )
    print('mlg',OSPL_dBA_mlg , 'dBA' )
    print('nlg',OSPL_dBA_nlg , 'dBA' )
    print('m strut',OSPL_dBA_mlg_strut, 'dBA'  )
    print('n strut',OSPL_dBA_nlg_strut, 'dBA' )

    print('overall dBA',OSPL_dBA_tot, 'dBA' )
    
    return  OSPL_dBA_tot #OSPL_dBA_flap, OSPL_dBA_slat,OSPL_dBA_wing,OSPL_dBA_mlg , OSPL_dBA_nlg ,OSPL_dBA_mlg_strut,OSPL_dBA_nlg_strut,
          
