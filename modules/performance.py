import numpy as np
from modules.fuel_consumption import fuel_flow

# requirement; take-off run 2000 meters @ sea level

# """
# Constants
# """
# rho = 1.225
# g = 9.81
# h_screen = 15.24
#
# """
# Variables
# """
# MTOW = 49391.28866 # kg
# thrust_takeoff_one_engine = 108.54E3 # N
# thrust_climb_out_one_engine = 108.54E3 * .85
# C_D = 0.04
# C_L_max = 1.4
# S = 113.8873926


def climb_gradient(T, D, m, g):
    y = (T-D)/m/g
    return y


def rate_of_climb(T, D, V, m, g):
    # from Lan, C.E. and Roskam, J., Airplane Aerodynamics & Performance, page 384
    RC = climb_gradient(T, D, m, g) * V
    return RC


def stall_speed(W, T, alpha_C_L_max, rho, C_L_max, S):
    theta_thrust = 0 # thrustline inclination
    V_stall = np.sqrt(2*(W-T*np.sin(alpha_C_L_max + theta_thrust))/(rho*C_L_max*S))
    return V_stall


def friction_coefficient(pressure_nose, weight, x_m, x_n, x_cg, h_cg):
    mu = ((pressure_nose / weight - 1) * (x_m - x_cg) + pressure_nose / weight * (x_cg - x_n))\
         / (h_cg * (1 - pressure_nose / weight))
    return mu


def take_off_field_length(rho, g, h_screen, MTOW, thrust_takeoff_one_engine, thrust_climb_out_one_engine, C_D, C_L_TO,
                          S, mu_TO): #MTOW in kg
    """
    Ground run
    """
    V_min = np.sqrt(MTOW * g / .5 / rho / S / C_L_TO)
    V_liftoff = 1.05 * V_min
    V_average = V_liftoff / np.sqrt(2)
    drag_ground_run_air = 0.5 * rho * V_average**2 * C_D * S
    average_lift = 0.5 * rho * V_average**2 * C_L_TO * S
    drag_ground_run_friction = mu_TO * (MTOW * g - average_lift)
    acceleration = (thrust_takeoff_one_engine - drag_ground_run_air - drag_ground_run_friction) / MTOW

    distance_ground = V_liftoff**2 / 2 / acceleration

    """
    Transition distance
    """
    drag_air = 0.5 * rho * V_liftoff ** 2 * C_D * S
    climb_angle = climb_gradient(thrust_takeoff_one_engine, drag_air, MTOW, g)
    # yamma = 0.9 * 0.295 - 0.3 / np.sqrt(AR)
    # same answer but we take climb_angle instead of yamma because why approximate if you can find it more exactly?

    x_transition = V_liftoff**2/0.15/g*np.sin(climb_angle)
    h_transition = V_liftoff**2/0.15/g*(1-np.cos(climb_angle))

    """
    Climb out distance
    """
    x_climb = 0
    if h_transition < h_screen:
        climb_out_angle = climb_gradient(thrust_climb_out_one_engine, drag_air, MTOW, g)
        x_climb = (h_screen - h_transition) / np.tan(climb_out_angle)
    x_total_airborne = x_transition + x_climb
    x_total = distance_ground + x_total_airborne
    # print(x_total)
    return x_total


# x_total = take_off_field_length(rho, g, h_screen, MTOW, thrust_takeoff_one_engine, thrust_climb_out_one_engine, C_D,
#                                 C_L_max, S, 0.05)
# print('x_total', x_total)


def landing_field_length(maximum_thrust, MTOW, g, h_screen, rho, S, C_L_LA, C_D, mu_LA):
    minimum_time_before_landing = 10*60 # s
    m_landing = MTOW - fuel_flow(0.7 * maximum_thrust) * minimum_time_before_landing
    V_min = np.sqrt(m_landing * g / .5 / rho / S / C_L_LA)
    V_approach = 1.3 * V_min
    delta_n = 0.10*g
    gamma_approach = 3/180*np.pi

    """
    Airborne distance
    """
    R = 1.3**2 * V_approach**2 / delta_n / g
    x_total_airborne = R * np.sin(gamma_approach) + (h_screen - (1 - np.cos(gamma_approach)) * R) / np.tan(gamma_approach)

    """
    Rotation distance
    """
    # assumption how long it takes to rotate;
    t_tr = 2.
    x_tr = t_tr * V_approach

    """
    Braking distance
    """
    reverse_thrust = 0.5 * maximum_thrust
    V_average = V_approach / np.sqrt(2)
    average_drag_air = 0.5 * rho * V_average ** 2 * C_D * S
    average_lift = 0.5 * rho * V_average ** 2 * C_L_LA * S
    drag_braking_friction = mu_LA * (m_landing * g - average_lift)
    x_brake = (m_landing*g)**2/g/S/rho*1.3**2/C_L_LA/(reverse_thrust+average_drag_air+drag_braking_friction)

    x_total = x_total_airborne + x_tr + x_brake
    return x_total
