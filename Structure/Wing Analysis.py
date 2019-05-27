# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:50:25 2019

@author: thong

Wing Analysis using Concept 1 values
Analyse canard
Analyse flap extension
Analyse wing tip extension
"""

"""
Parameters
Configuration 1
"""

b1 = 38.593274          # Span [m]
S1 = 135.403711         # Surface Area [m^2]
Cr1 = 5.359098          # Root Chord [m]
Ct1 = 1.657860          # Tip Chord [m]
MTOM1 = 58722.6         # lightest Maximum takeoff Mass [kg]

AR1 = b1**2/S1
TR1 = Ct1/Cr1

MTOW1 = MTOM1*9.81      # [N]
K = MTOW1/S1

"""
Parameters
Configuration 3
"""
MTOM2 = 68264.27        # Heaviest Maximum takeoff mass [kg]

MTOW2 = MTOM2*9.81      # [N]
S2 = MTOW2/K            # [m^2]

delta_S = S2-S1         # [m^2]

"""
Concept 1 (Canard)
"""
S_canard = delta_S
b_canard = (S_canard*AR1)**0.5
TR_canard = TR1
Cr_canard = 2*S_canard/((1+TR_canard)*b_canard)
Ct_canard = Cr_canard*TR_canard

"""
Concept 2
"""
S_concept2 = S2
b_concept2 = b1
AR_concept2 = b_concept2**2/S2
Ct_concept2 = Ct1
Cr_c


