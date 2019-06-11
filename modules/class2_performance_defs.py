import numpy as np
from Structure.Wing.isa import isa


def get_climb_gradient(thrust, drag, mass, g):
    y = (thrust - drag) / mass / g
    return y


def get_take_off_field_length(engine_failure, rho, g, h_screen, mass, thrust_one_engine, thrust_transition_setting, thrust_climb_out_setting, C_L, C_D, S, mu_TO, V_1_guess, V_1, reverse_thrust_factor):

    def get_acceleration(V_min, V_max, thrust_total):
        V_average = (V_max - V_min) / np.sqrt(2) + V_min
        average_lift = 0.5 * rho * V_average ** 2 * C_L * S
        drag_ground_run_air = 0.5 * rho * V_average ** 2 * C_D * S
        drag_ground_run_friction = mu_TO * (mass * g - average_lift)
        a = (thrust_total - drag_ground_run_air - drag_ground_run_friction) / mass
        return a

    def get_decision_speed(velocity):

        def give_try(velocity):
            acceleration_0 = get_acceleration(0, velocity, 2 * thrust_one_engine)
            acceleration_1 = get_acceleration(velocity, V_liftoff, thrust_one_engine)
            distance_try = velocity ** 2 / 2. / acceleration_0 + (V_liftoff ** 2 / 2. - velocity ** 2 / 2.) / acceleration_1

            if engine_failure:
                reverse_thrust = reverse_thrust_factor * thrust_one_engine
            else:
                reverse_thrust = reverse_thrust_factor * 2 * thrust_one_engine
            V_average = velocity / np.sqrt(2)
            average_drag_air = 0.5 * rho * V_average ** 2 * C_D * S
            average_lift = 0.5 * rho * V_average ** 2 * C_L * S
            drag_braking_friction = mu_TO * (mass * g - average_lift)
            a_brake = -1 / mass * (reverse_thrust + average_drag_air + drag_braking_friction)
            x_brake = -velocity ** 2 / 2 / a_brake

            distance_try+=x_brake

            return(distance_try)

        difference = 100
        while abs(difference)>5:
            if difference > 0:
                velocity += .1
            else:
                velocity -= .1

            nominal_distance = get_take_off_field_length(True, rho, g, h_screen, mass, thrust_one_engine,
                                                         thrust_transition_setting, thrust_climb_out_setting, C_L, C_D,
                                                         S,
                                                         mu_TO, False, velocity, reverse_thrust_factor)[0]
            try_distance = give_try(velocity)
            difference = nominal_distance - try_distance
        return velocity

    V_min = get_V_stall(mass, g, rho, S, C_L)
    V_liftoff = 1.05 * V_min  #1.05 CS-25

    """
    Ground run
    """
    if engine_failure:
        if V_1_guess:
            V_1 = get_decision_speed(V_1)
            acceleration_0 = get_acceleration(0, V_1, 2 * thrust_one_engine)
            acceleration_1 = get_acceleration(V_1, V_liftoff, thrust_one_engine)
        else:
            acceleration_0 = get_acceleration(0, V_1, 2 * thrust_one_engine)
            acceleration_1 = get_acceleration(V_1, V_liftoff, thrust_one_engine)
        distance_ground = V_1**2 / 2. / acceleration_0 + (V_liftoff**2/2. - V_1**2/2.) / acceleration_1
    else:
        acceleration = get_acceleration(0, V_liftoff, 2 * thrust_one_engine)
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

    else:
        gamma_2 = np.arccos(1-h_screen/(V_liftoff**2/0.15/g))
        distance_total_airborne = V_liftoff**2/0.15/g*np.sin(gamma_2)

    distance_total = distance_ground + distance_total_airborne

    return distance_total, V_liftoff


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
    return x_total, V_approach


def get_V_stall(mass, g, rho, S, C_L):
    V_stall = np.sqrt(mass * g * 2 / rho / S / C_L)
    return V_stall


def get_fuel_consumption(thrust, distance, velocity):
    # only true for PW1525G
    # approximately linear relation between thrust and fuel consumption
    # from linear regression;
    fuel_flow = 7E-6 * thrust + 0.0147                                                             # kg/s
    # R^2 = 0.9981
    total_fuel_used = distance / velocity * fuel_flow
    return fuel_flow, total_fuel_used


def get_energy_height():
    V = np.linspace(0, 250, 100)
    h = np.linspace(0, 11277.6, 100)
    z = []
    for i in range(0, 100, 1):
        list = []
        for j in range(0, 100, 1):
            list.append(h[i]+V[j]**2/2/9.80655)
        z.append(list)
    return V,h,z


def get_rate_of_climb(thrust, wing_surface_area, drag_coefficient, mass, g):
    V = np.linspace(0, 250, 100)
    h = np.linspace(0, 11277.6, 100)
    z = []
    for i in range(0, 100, 1):
        list = []
        for j in range(0, 100, 1):
            list.append((thrust-get_thrust_required(isa(h[i])[2], V[j], wing_surface_area, drag_coefficient))*V[j]/mass/g)
        z.append(list)
    return V, h, z


def get_range_breguet(mass_initial, mass_final, velocity, cj, LD_cruise, g):
    range_cruise = velocity/cj/g*LD_cruise*np.log(mass_initial/mass_final)
    return range_cruise


def get_fuel_burned_breguet(mass_initial, range_cruise, velocity, cj, g, LD_cruise):
    mass_final = mass_initial / np.exp(range_cruise * cj * g / velocity / LD_cruise)
    mass_fuel = mass_initial - mass_final
    return mass_fuel


def get_thrust_required(rho, velocity, S_wing, C_D_cruise):
    T = 0.5 * rho * velocity**2 * S_wing * C_D_cruise     # N
    return T
