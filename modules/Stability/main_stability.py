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

