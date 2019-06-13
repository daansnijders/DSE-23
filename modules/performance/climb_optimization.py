# maximum ROC, because it uses much fuel, minimum horizontal distance or minimum fuel don't come close
# chopped up into parts
# from inputs.concept_1 import thrust_max, H_m, V_cruise, S, CDcruise, MTOW, g
from Structure.Wing.isa import isa
# from modules.class2_performance_defs import *
# import matplotlib.pyplot as plt


def get_climb_optimization():
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx


    """
    inputs
    """
    steps = 50
    take_off_velocity = 70.
    i = 0
    mass = MTOW
    engines_operative = 2
    thrust_setting = 0.9
    climb_thrust = thrust_max * thrust_setting
    CD = CDcruise

    """
    analysis
    """


    a, b, rate_of_climb = get_2d_rate_of_climb(climb_thrust, engines_operative, S, CD[i], MTOW[i], g)
    cp = plt.contour(a**2/2/9.80655, b, rate_of_climb, 2*steps+1)
    plt.clabel(cp, inline=False, fontsize=10)
    # inline = True looks prettier but destroys the program as gradient can not be found for inline pieces


    """"
    finding raaklijn
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
        if value <= rico+0.1:
            if value >= rico-0.1:
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
        z = get_rate_of_climb(engines_operative, thrust_max, h_optimal[a], np.sqrt(V_optimal[a]*2*g), S, CD[i], mass[i], g)
        z_optimal.append(z)
    # h_optimal.append(take_off_velocity**2/2/g)
    # fn_roots = np.roots(rate_of_climb_fn-take_off_velocity**2/2/g)
    # V_optimal.append(np.real((fn_roots[np.isreal(fn_roots)])[0]))
    # z_optimal.append(get_rate_of_climb(engines_operative, thrust_max, take_off_velocity**2/2/g, np.sqrt(V_optimal[-1]*2*g), S, CD[i], mass[i], g))

    h_optimal = list(np.flip(h_optimal))
    V_optimal = list(np.flip(V_optimal))
    z_optimal = list(np.flip(z_optimal))
    # h_optimal.append(H_m)
    # V_optimal.append(V_cruise**2/2/g)
    # z_optimal.append(get_rate_of_climb(engines_operative, climb_thrust, V_optimal[-1], np.sqrt(V_optimal[-1]*2*g), S, CD[i], mass[i], g))

    rate_of_climb_fn = np.poly1d(np.polyfit(V_optimal, h_optimal, 3))

    climb_time = 0
    distance = 0
    fuel_mass_climb = 0
    for a in range(1, len(h_optimal)):
        time_add = (h_optimal[a]-h_optimal[a-1])/((z_optimal[a]+z_optimal[a-1])/2)
        climb_time += time_add
        distance += (h_optimal[a]-h_optimal[a-1])/np.tan(np.arcsin(((z_optimal[a]+z_optimal[a-1])/2)/((V_optimal[a]+V_optimal[a-1])/2)))
        fuel_mass_climb += (h_optimal[a]-h_optimal[a-1])*(engines_operative*get_fuel_consumption(climb_thrust, 1, 1)[0])/((z_optimal[a]+z_optimal[a-1])/2)
    # plt.scatter(x_loc_list, y_loc_list)
    plt.scatter(V_optimal, h_optimal)
    polyfit_x = np.linspace(500, 2000, 1000)
    polyfit_y = rate_of_climb_fn(polyfit_x)
    plt.plot(polyfit_x, polyfit_y)
    x, y, energy_height = get_energy_height()
    plt.contour(x**2/2/9.80655, y, energy_height, steps+1)
    plt.xlim(left=0.0)
    plt.ylim(bottom=0.0, top=20000.)

    fuel_flow_climb = engines_operative*get_fuel_consumption(climb_thrust, 1, 1)[0]
    climb_final_velocity = V_optimal[-1]

    return fuel_flow_climb, fuel_mass_climb, climb_final_velocity
