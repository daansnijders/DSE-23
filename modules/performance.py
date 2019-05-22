import numpy as np

# requirement; take-off run 2000 meters @ sea level


def get_climb_gradient(T, D, m, g):
    y = (T-D)/m/g
    return y


def get_rate_of_climb( V, climb_gradient):
    # from Lan, C.E. and Roskam, J., Airplane Aerodynamics & Performance, page 384
    RC = climb_gradient * V
    return RC


def get_stall_speed(W, T, alpha_C_L_max, rho, C_L_max, S):
    theta_thrust = 0 # thrustline inclination
    V_stall = np.sqrt(2*(W-T*np.sin(alpha_C_L_max + theta_thrust))/(rho*C_L_max*S))
    return V_stall


def get_friction_coefficient(force_nose, mass, x_m, x_n, x_cg, h_cg, g):
    weight = mass * g
    mu = ((force_nose / weight - 1) * (x_m - x_cg) + force_nose / weight * (x_cg - x_n))\
         / (h_cg * (1 - force_nose / weight))
    return abs(mu)


def get_take_off_field_length(rho, g, h_screen, MTOW, thrust_takeoff_one_engine, thrust_climb_out_one_engine, C_D,
                              C_L_TO, S, mu_TO): #MTOW in kg
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
    # climb_angle = 8/180*np.pi
    climb_angle = get_climb_gradient(2*0.8*thrust_takeoff_one_engine, drag_air, MTOW, g)
    #
    # yamma = 0.9 * 0.295 - 0.3 / np.sqrt(AR)
    # same answer but we take climb_angle instead of yamma because why approximate if you can find it more exactly?

    x_transition = V_liftoff**2/0.15/g*np.sin(climb_angle)
    h_transition = V_liftoff**2/0.15/g*(1-np.cos(climb_angle))

    """
    Climb out distance
    """
    x_total_airborne = x_transition
    if h_transition < h_screen:
        climb_out_angle = get_climb_gradient(thrust_climb_out_one_engine, drag_air, MTOW, g)
        x_climb = (h_screen - h_transition) / np.tan(climb_out_angle)
        x_total_airborne = x_transition + x_climb
    x_total = distance_ground + x_total_airborne
    # print(x_total)
    return x_total


def get_m_landing(MTOW, maximum_thrust):
    minimum_time_before_landing = 20 * 60  # s
    m_landing = MTOW - get_fuel_flow(0.7 * maximum_thrust) * minimum_time_before_landing
    return m_landing


def get_landing_field_length(maximum_thrust, m_landing, g, h_screen, rho, S, C_L_LA, C_D, mu_LA):
    V_min = np.sqrt(m_landing * g / .5 / rho / S / C_L_LA)
    V_approach = 1.3 * V_min
    delta_n = 0.10*g
    gamma_approach = 3./180*np.pi

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
    a_brake = -1/m_landing*(reverse_thrust+average_drag_air+drag_braking_friction)
    x_brake = -V_approach**2/2/a_brake
    x_total = x_total_airborne + x_tr + x_brake
    return x_total


def get_fuel_flow(thrust):
    # only true for PW1525G
    # approximately linear relation between thrust and fuel consumption
    # from linear regression;
    ff = 7E-6 * thrust + 0.0147                                                             # kg/s
    # R^2 = 0.9981
    return ff


def get_cruise_thrust(rho, cruise_velocity, S_wing, C_D_cruise):
    T = 0.5 * rho * cruise_velocity**2 * S_wing * C_D_cruise     # N
    return T


def get_cruise_fuel(thrust_cruise, cruise_range, cruise_velocity):
    fuel_flow_cruise = get_fuel_flow(thrust_cruise)
    total_fuel_used_cruise = cruise_range / cruise_velocity * fuel_flow_cruise
    return fuel_flow_cruise, total_fuel_used_cruise


def get_V_min(W, g, rho, S, C_L):
    V_min = np.sqrt(W * g / .5 / rho / S / C_L)
    return V_min


def specific_fuel_consumption(thrust, fuel_flow):
    poundforce = thrust * 0.224809
    fuel_flow_pound_per_hour = fuel_flow * 2.20462 / 3600
    fuel_flow_pound_per_second = fuel_flow * 2.20462
    cj_si = fuel_flow / thrust
    cj_lbhour = fuel_flow_pound_per_hour/poundforce
    cj_lbsec = fuel_flow_pound_per_second/poundforce
    return cj_si, cj_lbhour, cj_lbsec

