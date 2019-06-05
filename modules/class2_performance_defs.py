import numpy as np


def get_climb_gradient(thrust, drag, mass, g):
    y = (thrust - drag) / mass / g
    return y


def get_take_off_field_length(engine_failure, rho, g, h_screen, mass, thrust_one_engine, thrust_transition_setting, thrust_climb_out_setting, C_L, C_D, S, mu_TO):

    def get_acceleration(V_max, thrust_total):
        V_average = V_max / np.sqrt(2)
        average_lift = 0.5 * rho * V_average ** 2 * C_L * S
        drag_ground_run_air = 0.5 * rho * V_average ** 2 * C_D * S
        drag_ground_run_friction = mu_TO * (mass * g - average_lift)
        a = (thrust_total - drag_ground_run_air - drag_ground_run_friction) / mass
        return a

    V_min = get_V_stall(mass, g, rho, S, C_L)
    V_liftoff = 1.05 * V_min

    """
    Ground run
    """
    if engine_failure:
        V_1 = 40.                                                      # todo; determine V1
        acceleration_0 = get_acceleration(V_1, 2 * thrust_one_engine)
        acceleration_1 = get_acceleration(V_liftoff, thrust_one_engine)
        distance_ground = V_1**2 / 2. / acceleration_0 + (V_liftoff**2/2. - V_1**2/2.) / acceleration_1
    else:
        acceleration = get_acceleration(V_liftoff, 2 * thrust_one_engine)
        distance_ground = V_liftoff**2 / 2. / acceleration
    """
    Transition distance
    """
    drag_air = 0.5 * rho * V_liftoff ** 2 * C_D * S
    # climb_gradient = 8/180*np.pi
    if engine_failure:  # todo; review climb gradient
        climb_gradient = get_climb_gradient(thrust_transition_setting * thrust_one_engine, drag_air, mass, g)
    else:
        climb_gradient = get_climb_gradient(2 * thrust_transition_setting * thrust_one_engine, drag_air, mass, g)
    distance_transition = V_liftoff**2/0.15/g*np.sin(climb_gradient)
    h_transition = V_liftoff**2/0.15/g*(1-np.cos(climb_gradient))

    """
    Climb out distance
    """
    distance_total_airborne = distance_transition
    if h_transition < h_screen:
        if engine_failure:
            climb_out_gradient = get_climb_gradient(thrust_climb_out_setting * thrust_one_engine, drag_air, mass, g)
        else:
            climb_out_gradient = get_climb_gradient(2 * thrust_climb_out_setting * thrust_one_engine, drag_air, mass, g)
        distance_climb = (h_screen - h_transition) / np.tan(climb_out_gradient)
        distance_total_airborne = distance_transition + distance_climb
    distance_total = distance_ground + distance_total_airborne
    return distance_total


def get_friction_coefficient(force_nose, mass, x_m, x_n, x_cg, h_cg, g):
    weight = mass * g
    mu = ((force_nose / weight - 1) * (x_m - x_cg) + force_nose / weight * (x_cg - x_n))\
         / (h_cg * (1 - force_nose / weight))
    return abs(mu)


def get_landing_field_length(engine_failure, rho, g, h_screen, mass, thrust_one_engine, C_L, C_D, S, mu_LA, reverse_thrust_factor):
    V_min = np.sqrt(mass * g / .5 / rho / S / C_L)
    V_approach = 1.23 * V_min
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
    if engine_failure:
        reverse_thrust = reverse_thrust_factor * thrust_one_engine
    else:
        reverse_thrust = reverse_thrust_factor* 2 * thrust_one_engine
    V_average = V_approach / np.sqrt(2)
    average_drag_air = 0.5 * rho * V_average ** 2 * C_D * S
    average_lift = 0.5 * rho * V_average ** 2 * C_L * S
    drag_braking_friction = mu_LA * (mass * g - average_lift)
    a_brake = -1 / mass * (reverse_thrust + average_drag_air + drag_braking_friction)
    x_brake = -V_approach**2/2/a_brake
    x_total = x_total_airborne + x_tr + x_brake
    return x_total


def get_V_stall(mass, g, rho, S, C_L):
    V_stall = np.sqrt(mass * g * 2 / rho / S / C_L)
    return V_stall