# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:46:51 2019

@author: Lisa
"""
import numpy as np
from math import *
from modules.sustainability.noise_defs import *


req_approach_EPNdB=89.9
req_lateral_TO_EPNdB=86
experimental_dBA_to_EPNdB=1.15503 #false


limit_epndb=[get_the_limit_for_lateral(MTOW[i]) for i in range(3)]










number_of_bands=43
bandnumbers=list(range(1,number_of_bands+1))
centrefreq=[10**(bandnumbers[i]/10) for i in range(len(bandnumbers))]
lowerfreq=[2**(-1/6)*centrefreq[i] for i in range(len(bandnumbers))]
upperfreq=[2**(1/6)*centrefreq[i] for i in range(len(bandnumbers))]
freq_delta=[upperfreq[i]-lowerfreq[i] for i in range(len(bandnumbers))]

