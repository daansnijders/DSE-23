from inputs.concept_1 import ft_to_m, thrust_max

"""
inputs
"""
i = 2   # configuration selection
h_screen_to = 35 * ft_to_m                                  # [m]
h_screen_la = 50 * ft_to_m
reverse_thrust_factor = 0.45
engine_failure = False
cj = 0.790/thrust_max #kg/s/N