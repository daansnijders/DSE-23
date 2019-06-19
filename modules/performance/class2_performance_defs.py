import numpy as np
from Structure.Wing.isa import isa
import matplotlib.pyplot as plt


def get_climb_gradient(thrust, drag, mass, g):
    y = (thrust - drag) / mass / g
    return y


def get_take_off_field_length(engine_failure, rho, g, h_screen, mass, thrust_one_engine, thrust_transition_setting, thrust_climb_out_setting, C_L, C_D, S, mu_TO, V_1_guess, V_1, reverse_thrust_factor, mu_LA):

    def get_acceleration(V_min, V_max, thrust_total):
        V_average = (V_max - V_min) / np.sqrt(2) + V_min
        average_lift = 0.5 * rho * V_average ** 2 * C_L * S
        drag_ground_run_air = 0.5 * rho * V_average ** 2 * C_D * S
        drag_ground_run_friction = mu_TO * (mass * g - average_lift)
        a = (thrust_total - drag_ground_run_air - drag_ground_run_friction) / mass
        return a

    def get_decision_speed(velocity):

        def give_try(decision_speed):
            acceleration = get_acceleration(0, decision_speed, 2 * thrust_one_engine)
            accelerate_stop_distance = decision_speed ** 2 / 2. / acceleration

            reverse_thrust = reverse_thrust_factor * thrust_one_engine
            average_velocity = decision_speed / np.sqrt(2)
            average_drag_air = 0.5 * rho * average_velocity ** 2 * C_D * S
            average_lift = 0.5 * rho * average_velocity ** 2 * C_L * S
            drag_braking_friction = mu_LA * (mass * g - average_lift)
            a_brake = -1 / mass * (reverse_thrust + average_drag_air + drag_braking_friction)
            x_brake = -decision_speed ** 2 / 2 / a_brake
            accelerate_stop_distance += x_brake
            return accelerate_stop_distance

        difference = 100
        count = 0
        # nominal_list = []
        # try_list = []
        # difference_list = []
        # count_list = []
        while abs(difference)>5:
            if abs(difference)>60:
                step_size = 1
            elif abs(difference)>40:
                step_size = 0.5
            else:
                step_size = 0.1
            count += 1
            if difference > 0:
                velocity += step_size

            else:
                velocity -= step_size
            nominal_distance = get_take_off_field_length(True, rho, g, h_screen, mass, thrust_one_engine,
                                                         thrust_transition_setting, thrust_climb_out_setting, C_L, C_D,
                                                         S, mu_TO, False, velocity, reverse_thrust_factor, mu_LA)[0]
            try_distance = give_try(velocity)
            difference = nominal_distance - try_distance
            if count>1000:
                print('ALERT - decision speed does not converge, check input values')
                velocity = 30
                break
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
        if V_1_guess:
            V_1 = get_decision_speed(V_1)
            acceleration_0 = get_acceleration(0, V_liftoff, 2 * thrust_one_engine)
            acceleration_1 = get_acceleration(V_1, V_liftoff, 2 * thrust_one_engine)
        else:
            acceleration_0 = get_acceleration(0, V_1, 2 * thrust_one_engine)
            acceleration_1 = get_acceleration(V_1, V_liftoff, 2 * thrust_one_engine)
        distance_ground = V_1 ** 2 / 2. / acceleration_0 + (V_liftoff ** 2 / 2. - V_1 ** 2 / 2.) / acceleration_1
        # acceleration = get_acceleration(0, V_liftoff, 2 * thrust_one_engine)
        # distance_ground = V_liftoff**2 / 2. / acceleration
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
    if h_transition > h_screen:
        distance_total_airborne = distance_transition
    else:
        if engine_failure:
            climb_out_gradient = get_climb_gradient(thrust_climb_out_setting * thrust_one_engine, drag_air, mass, g)
        else:
            climb_out_gradient = get_climb_gradient(2 * thrust_climb_out_setting * thrust_one_engine, drag_air, mass, g)
        distance_climb = (h_screen - h_transition) / np.tan(climb_out_gradient)
        distance_total_airborne = distance_transition + distance_climb
    # if h_transition < h_screen:
    #     if engine_failure:
    #         climb_out_gradient = get_climb_gradient(thrust_climb_out_setting * thrust_one_engine, drag_air, mass, g)
    #     else:
    #         climb_out_gradient = get_climb_gradient(2 * thrust_climb_out_setting * thrust_one_engine, drag_air, mass, g)
    #     distance_climb = (h_screen - h_transition) / np.tan(climb_out_gradient)
    #     distance_total_airborne = distance_transition + distance_climb
    #
    # else:
    #     gamma_2 = np.arccos(1-h_screen/(V_liftoff**2/0.15/g))
    #     distance_total_airborne = V_liftoff**2/0.15/g*np.sin(gamma_2)

    distance_total = distance_ground + distance_total_airborne
    return distance_total, V_liftoff, V_1


def get_take_off_field_length_alt(altitude, CD_TO, CL_TO, friction_coefficient_to, g, mass, S, screen_height_to, thrust_max):
    air_density = isa(altitude)[2]
    stall_velocity = get_V_stall(mass, g, air_density, S, CL_TO)
    take_off_velocity = 1.05 * stall_velocity
    average_velocity = take_off_velocity / np.sqrt(2)
    average_acceleration = (thrust_max * 2 - friction_coefficient_to * (mass * g - CL_TO * 1 / 2 * air_density * average_velocity ** 2 * S) - CD_TO * 1 / 2 * air_density * average_velocity ** 2 * S) / mass
    distance_ground = take_off_velocity ** 2 / 2 / average_acceleration

    # confirmed it works semi, deltas = 40

    radius = take_off_velocity ** 2 / g / 0.15
    flight_path_angle = (2 * thrust_max - 1 / 2 * air_density * take_off_velocity ** 2 * S * CD_TO) / mass / g
    distance_transition = radius * flight_path_angle
    height_transition = 1 / 2 * distance_transition * flight_path_angle

    distance_climb = (screen_height_to - height_transition) / np.tan(flight_path_angle)
    if height_transition > screen_height_to:
        distance = distance_transition + distance_ground
    else:
        distance = distance_transition + distance_ground + distance_climb
    return distance, take_off_velocity, 40


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
    height_intervals = 10
    V = np.linspace(0, 300, height_intervals)
    h = np.linspace(0, 20000, height_intervals)
    z = []
    for i in range(0, height_intervals, 1):
        list = []
        for j in range(0, height_intervals, 1):
            list.append(h[i]+V[j]**2/2/9.80655)
        z.append(list)
    return V,h,z


def get_2d_rate_of_climb(thrust_max, engines_operative, wing_surface_area, drag_coefficient, mass, g):
    steps_V = 100
    steps_h = 1000
    V = np.linspace(0, 300, steps_V)
    h = np.linspace(0, 20000, steps_h)
    z = []
    for i in range(0, steps_h, 1):
        list = []
        for j in range(0, steps_V, 1):
            list.append(get_rate_of_climb(engines_operative, thrust_max, h[i], V[j], wing_surface_area, drag_coefficient, mass, g))
        z.append(list)

    return V, h, z


def get_rate_of_climb(engines_operative, thrust_max, h, V, wing_surface_area, drag_coefficient, mass, g):
    thrust_available = engines_operative * get_thrust_available(thrust_max, h, V)
    z = ((thrust_available - get_thrust_required(isa(h)[2], V, wing_surface_area, drag_coefficient)) * V / mass / g)
    return z


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


def get_thrust_available(thrust_max, altitude, velocity):
    # assume linear relation with density
    n = 1.0
    altitude_effect = (isa(altitude)[2] / isa(0)[2])**n

    pressure_ratio = 38.7
    fuel_flow = 0.790
    p0 = 101325


    velocity_effect = 1

    thrust = thrust_max * altitude_effect * velocity_effect
    return thrust


def get_climb_optimization(mass_climb_initial, thrust_max, CD_climb, S, g, H_m, V_cruise, thrust_setting):
    import matplotlib.pyplot as plt
    import numpy as np

    np.seterr(divide='ignore', invalid='ignore')

    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx

    """
    inputs
    """
    steps = 250
    mass = mass_climb_initial
    engines_operative = 2
    climb_thrust = thrust_max * thrust_setting

    """
    analysis
    """
    plt.figure()
    a, b, rate_of_climb = get_2d_rate_of_climb(climb_thrust, engines_operative, S, CD_climb, mass_climb_initial, g)
    cp = plt.contour(a**2/2/9.80655, b, rate_of_climb, 2*steps+1)
    plt.clabel(cp, inline=False, fontsize=10)
    # inline = True looks prettier but destroys the program as gradient can not be found for inline pieces

    """"
    finding tangent line
    """
    x_loc_list = []
    y_loc_list = []
    for a in range(1, len(cp.collections)-1):
        collections = cp.collections[a]
        x_list = []
        y_list = []
        for b in range(0, len(collections.get_paths())):
            paths = collections.get_paths()[b]
            v = paths.vertices
            x = v[:,0]
            y = v[:,1]
            x_list.extend(x)
            y_list.extend(y)
        gradient = np.gradient(y_list, x_list)
        rico = -1.
        value, index = find_nearest(gradient, rico)
        if rico-0.1 < value < rico+0.1:
            x_loc = x_list[index]
            y_loc = y_list[index]
            x_loc_list.append(x_loc)
            y_loc_list.append(y_loc)

    # filtering out points below cruise altitude
    below_cruise_altitude = [item <= H_m for item in y_loc_list]
    h_optimal = [d for (d, remove) in zip(y_loc_list, below_cruise_altitude) if remove]
    V_optimal = [d for (d, remove) in zip(x_loc_list, below_cruise_altitude) if remove]
    z_optimal = []
    for a in range(0, len(h_optimal)):
        z = get_rate_of_climb(engines_operative, thrust_max, h_optimal[a], np.sqrt(V_optimal[a]*2*g), S, CD_climb, mass, g)
        z_optimal.append(z)
    # h_optimal.append(take_off_velocity**2/2/g)
    # fn_roots = np.roots(rate_of_climb_fn-take_off_velocity**2/2/g)
    # V_optimal.append(np.real((fn_roots[np.isreal(fn_roots)])[0]))
    # z_optimal.append(get_rate_of_climb(engines_operative, thrust_max, take_off_velocity**2/2/g, np.sqrt(V_optimal[-1]*2*g), S, CD[i], mass[i], g))

    h_nox = 3000 * 0.3048

    h_optimal = list(np.flip(h_optimal, 0))
    V_optimal = list(np.flip(V_optimal, 0))
    z_optimal = list(np.flip(z_optimal, 0))


    # rate_of_climb_fn = np.poly1d(np.polyfit(V_optimal, h_optimal, 2))

    # z_optimal.append(get_rate_of_climb(engines_operative, climb_thrust, V_optimal[-1], np.sqrt(V_optimal[-1] * 2 * g), S, CD_climb[i], mass[i], g))
    # h_optimal.append(H_m)
    # V_optimal.append(V_cruise**2/2/g)

    climb_time = 0
    distance = 0
    fuel_mass_climb = 0
    fuel_mass_3000 = 0
    for a in range(1, len(h_optimal)):
        time_add = ((h_optimal[a]+V_optimal[a]**2/2/g)-(h_optimal[a-1]+V_optimal[a-1]**2/2/g))/((z_optimal[a]+z_optimal[a-1])/2)
        climb_time += time_add
        distance += ((h_optimal[a]+V_optimal[a]**2/2/g)-(h_optimal[a-1]+V_optimal[a-1]**2/2/g))/np.tan(np.arcsin(((z_optimal[a]+z_optimal[a-1])/2)/((V_optimal[a]+V_optimal[a-1])/2)))
        fuel_add = ((h_optimal[a]+V_optimal[a]**2/2/g)-(h_optimal[a-1]+V_optimal[a-1]**2/2/g))*(engines_operative*get_fuel_consumption(climb_thrust, 1, 1)[0])/((z_optimal[a]+z_optimal[a-1])/2)
        fuel_mass_climb += fuel_add
        if h_optimal[a] < h_nox:
            fuel_mass_3000 += fuel_add
    # plt.scatter(x_loc_list, y_loc_list)
    plt.scatter(V_optimal, h_optimal)
    # polyfit_x = np.linspace(500, 2000, 1000)
    # polyfit_y = rate_of_climb_fn(polyfit_x)
    # plt.plot(polyfit_x, polyfit_y)
    x, y, energy_height = get_energy_height()
    plt.contour(x**2/2/9.80655, y, energy_height, steps+1)
    plt.xlim(left=0.0)
    plt.ylim(bottom=0.0, top=20000.)

    fuel_flow_climb = engines_operative*get_fuel_consumption(climb_thrust, 1, 1)[0]
    climb_final_velocity = V_optimal[-1]

    return fuel_flow_climb, fuel_mass_climb, climb_final_velocity, distance, fuel_mass_3000


def get_descent(max_altitude, cruise_velocity, approach_velocity, engines_operative, thrust_descent, S, drag_coefficient, mass, g):
    fuel_flow_descent = get_fuel_consumption(thrust_descent, 1, 1)[0]

    pieces = 50
    descent_time = 0
    descent_distance = 0
    fuel_mass_descent = 0
    ROC_list = []
    altitude_list = np.linspace(max_altitude, 1, pieces)
    velocity_list = np.linspace(cruise_velocity, approach_velocity, pieces)
    for index in range(0, len(altitude_list)):
        altitude = altitude_list[index]
        velocity = velocity_list[index]
        ROC = get_rate_of_climb(engines_operative, thrust_descent, altitude, velocity, S, drag_coefficient, mass, g)
        ROC_list.append(ROC)
    for a in range(1, len(altitude_list)):
        if ROC_list[a] < 0:
            time_add = ((altitude_list[a]+velocity_list[a]**2/2/g)-(altitude_list[a-1]+velocity_list[a-1]**2/2/g))/((ROC_list[a]+ROC_list[a-1])/2)
            descent_time += time_add
            descent_distance += ((altitude_list[a]+velocity_list[a]**2/2/g)-(altitude_list[a-1]+velocity_list[a-1]**2/2/g))/np.tan(np.arcsin(((ROC_list[a]+ROC_list[a-1])/2)/((velocity_list[a]+velocity_list[a-1])/2)))
        else:
            break
        fuel_mass_descent = descent_time * (engines_operative*fuel_flow_descent)
    return fuel_mass_descent, descent_distance
