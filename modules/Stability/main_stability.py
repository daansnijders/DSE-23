# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:56:29 2019

@author: Stijn
"""
import numpy as np
import matplotlib.pyplot as plt

from inputs.concept_1 import *
from inputs.constants import *
from modules.main_class2 import *

from modules.Stability.canard import *

from modules.Stability.empennage import *
from modules.Stability.cg_weight_loadingdiagram import *
from modules.Stability.cg_weight_config1 import *


e2 = empennage(1, (11.78+0.25*3.8), 3.82, 4.90, 0.3835, 21.72, 16., 93.5, 3.8, 1., 11.78, -0.3, 1.6, x_cg_max, -0.5838, )



#c2 = canard(weight_pass,2, 1.3)        
#c3 = canard(weight_pass,3, 1.3)
#
#print('Canard: ' + str(c3.F_c / c2.F_c *100 - 100))
#print('Main w: ' + str(c3.F_w / c2.F_w *100 - 100))
#print('H tail: ' + str(c3.F_h / c2.F_h *100 - 100))