import numpy as np
# requirement; take-off run 2000 meters @ sea level

"""
Constants
"""
rho = 1.225
g = 9.81
h_screen = 10.

"""
Variables
"""
MTOW = 49391.28866 # kg
thrust_takeoff = 108.54E3 # N
C_D = 0.02
C_L_max = 0.7813
S = 113.8873926
AR = 11.

"""
Ground run
"""
V_liftoff = np.sqrt(MTOW * g / .5 / rho / S / C_L_max)
V_average = V_liftoff / np.sqrt(2)
drag_ground_run_air = 0.5 * rho * V_average** 2 * C_D * S
drag_ground_run_friction = 0

acceleration = (thrust_takeoff - drag_ground_run_air - drag_ground_run_friction) / MTOW

distance_ground = V_liftoff**2 / 2 / acceleration

"""
Transition distance
"""
drag_air = 0.5 * rho * V_liftoff ** 2 * C_D * S
climb_angle = (thrust_takeoff - drag_air)/MTOW/9.81
yamma = 0.9 * 0.295 - 0.3 / np.sqrt(AR)
# same answer but we take climb_angle instead of yamma because why approximate if you can find it more exactly?

x_transition = V_liftoff**2/0.15/g*np.sin(climb_angle)
h_transition = V_liftoff**2/0.15/g*(1-np.cos(climb_angle))

"""
Climb out distance
"""
x_climb = (h_screen - h_transition) / np.tan(climb_angle)

x_total_airborne = x_transition - x_climb
x_total = distance_ground + x_total_airborne
print(x_total)
