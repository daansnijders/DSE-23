import matplotlib.pyplot as plt
import pandas as pd

from Structure.Wing.isa import isa
from inputs.concept_1 import ft_to_m, OEW, MTOW, thrust_max, S, g, M_payload, M_fuel, A, e, CD0, R, V_cruise, LoverD, CDcruise
from modules.class2_performance_defs import *
from inputs.constants import H_m

"""
inputs
"""
i = 0   # configuration selection
h_screen_to = 35 * ft_to_m                                  # [m]
h_screen_la = 50 * ft_to_m
reverse_thrust_factor = 0.45
engine_failure = False
cj = 0.790/thrust_max #kg/s/N
cj_retard = cj / 2.832545035E-5
thrust_transition_setting = 1.
thrust_climb_out_setting = 1.  # todo; reconsider this value (sometimes 0.85)
C_L_to = 1.9  # todo; lift and drag of canard (update C_L_TO, C_D_TO, S)
C_D_to = CD0[i] + 0.015 + 0.02 + C_L_to ** 2 / np.pi / A / e  # todo; lift and drag of canard (update C_L_TO, C_D_TO, S, C_D lift induced)
correction_factor_to = 1.15  # from CS-25
# friction_coefficient = get_friction_coefficient(P_nw[i], mass, x_mlg[i], x_nlg[i], x_cg[i], z_cg - z_mlg[i], g)
friction_coefficient_to = 0.05  # from Raymer table 17.1

C_L_la = 2.4                    # todo; lift and drag of canard (update C_L_LA, C_D_LA, S)
C_D_la = CD0[i] + 0.015 + 0.02 + C_L_la**2/np.pi/A/e                  # todo; lift and drag of canard (update C_L_LA, C_D_LA, S)
# friction_coefficient = get_friction_coefficient(P_nw[i], mass, x_mlg[i], x_nlg[i], x_cg[i], z_cg - z_mlg[i], g)
friction_coefficient_la = 0.45      # from Raymer table 17.1
correction_factor_la = 1.          # todo; find in CS-25??? look it up for landing


"""
take-off performance
"""


def analyze_take_off_performance():
    # todo; determine decision speed (need brake design)

    max_airport_altitude = 2000.                            # [m]
    altitude_resolution = 5
    airport_altitude_list = np.linspace(0, max_airport_altitude, altitude_resolution)

    take_off_field_length = []

    mass_resolution = 50
    mass_list = np.linspace(OEW[i], MTOW[i], mass_resolution)   # todo; incorrect

    for altitude in airport_altitude_list:
        density = isa(altitude)[2]
        take_off_field_length_list = []
        take_off_velocity_list = []
        for mass in mass_list:
            length, take_off_velocity = get_take_off_field_length(engine_failure, density, g, h_screen_to, mass,
                                                                  thrust_max, thrust_transition_setting,
                                                                  thrust_climb_out_setting, C_L_to, C_D_to, S,
                                                                  friction_coefficient_to, True, 50, reverse_thrust_factor)
            take_off_field_length_list.append(length * correction_factor_to)
            take_off_velocity_list.append(take_off_velocity)
        take_off_field_length.append([mass_list, take_off_field_length_list, take_off_velocity_list])

    # plotting
    plt.figure()
    for select in take_off_field_length:
        h = airport_altitude_list[take_off_field_length.index(select)]
        plt.plot(select[1], select[0], label='%a [m]' % h)

    plt.legend()
    engines_used = 2 - engine_failure

    plt.title('Take-off field length, %a engine(s) operative' % engines_used)

    plt.axhline(y=OEW[i], linestyle=':')
    plt.text(2000, OEW[i], 'OEW')
    plt.axhline(y=MTOW[i], linestyle=':')
    plt.text(2000, MTOW[i], 'MTOW')
    plt.axhline(y=OEW[i]+M_payload[i], linestyle=':')
    plt.text(2000, OEW[i] + M_payload[i], 'OEW+M_payload')
    plt.axhline(y=OEW[i]+M_fuel[i], linestyle=':')
    plt.text(2000, OEW[i]+M_fuel[i], 'OEW+M_fuel')

    plt.axvline(2000, linestyle=':')

    plt.ylabel('Mass [kg]')
    plt.xlabel('Take-off field length [m]')
    plt.close()
    return take_off_field_length, airport_altitude_list


# take_off_field_length = analyze_take_off_performance()[0]

"""
landing performance
"""


def analyze_landing_performance():
    # inputs
    max_airport_altitude = 1000.  # [m]
    altitude_resolution = 5
    airport_altitude_list = np.linspace(0, max_airport_altitude, altitude_resolution)

    landing_field_length = []

    mass_resolution = 50
    mass_list = np.linspace(OEW[i], MTOW[i], mass_resolution)   # todo; review this, it's incorrect, will not land with MTOW

    for altitude in airport_altitude_list:
        density = isa(altitude)[2]
        landing_field_length_list = []
        approach_velocity_list = []
        for mass in mass_list:
            x, V_approach = get_landing_field_length(engine_failure, density, g, h_screen_la, mass, thrust_max,C_L_la,
                                                     C_D_la, S, friction_coefficient_la, reverse_thrust_factor)
            landing_field_length_list.append(x * correction_factor_la)
            approach_velocity_list.append(V_approach)
        landing_field_length.append([mass_list, landing_field_length_list, approach_velocity_list])

    # plotting
    plt.figure()
    for select in landing_field_length:
        h = airport_altitude_list[landing_field_length.index(select)]
        plt.plot(select[1], select[0], label='%a [m]' % h)

    plt.legend()
    engines_used = 2 - engine_failure
    plt.title('Landing field length, %a engine(s) operative' % engines_used)
    # plt.plot([0, max(select[1])], [MTOW[i], MTOW[i]], color='C0', linestyle=':')
    # plt.plot([0, max(select[1])], [OEW[i], OEW[i]], color='C0', linestyle=':')
    plt.ylabel('Mass [kg]')
    plt.xlabel('Landing field length [m]')
    plt.close()
    return landing_field_length


landing_field_length = analyze_landing_performance()


"""
fuel economy
"""


def analyze_fuel_consumption(MTOW, take_off_velocity):
    fuel_consumption = pd.DataFrame(data=None, columns=['fuel_flow', 'fuel_mass'])
    pick_breguet = True
    # engine_failure = False
    engines_operative = 2 - engine_failure
    """
    engine startup/taxi
    """
    fuel_fraction_engine_startup = 0.990
    fuel_fraction_taxi = 0.990

    fuel_mass_engine_startup = (1-fuel_fraction_engine_startup)*MTOW[i]
    fuel_mass_taxi = (1-fuel_fraction_taxi)*MTOW[i]
    fuel_consumption.loc['engine_startup'] = [np.NaN, fuel_mass_engine_startup]
    fuel_consumption.loc['taxi'] = [np.NaN, fuel_mass_taxi]

    """
    take-off !!!!! FOR SEA LEVEL TAKE-OFF !!!!
    """
    # assume full thrust during entire take-off


    take_off_field_length, take_off_velocity = get_take_off_field_length(engine_failure, isa(0)[2], g, h_screen_to, MTOW[i], thrust_max, thrust_transition_setting, thrust_climb_out_setting, C_L_to, C_D_to, S, friction_coefficient_to, True, 80, reverse_thrust_factor)
    fuel_flow_take_off, fuel_mass_take_off = get_fuel_consumption(engines_operative*thrust_max, take_off_field_length, take_off_velocity/np.sqrt(2))
    fuel_consumption.loc['take_off'] = [fuel_flow_take_off, fuel_mass_take_off]

    mass = MTOW[i] - fuel_consumption.loc['take_off']['fuel_mass'] - fuel_consumption.loc['engine_startup']['fuel_mass'] - fuel_consumption.loc['taxi']['fuel_mass']

    """
    climb
    """

    fuel_fraction_climb = 0.980
    fuel_mass_climb = (1 - fuel_fraction_climb) * MTOW[i]
    fuel_consumption.loc['climb'] = ['variable', fuel_mass_climb]
    mass -= fuel_mass_climb


    """
    cruise breguet
    """
    fuel_mass_cruise_breguet = get_fuel_burned_breguet(mass, R[i], V_cruise, cj, g, LoverD[i])  # todo; update cj
    fuel_flow_cruise_breguet = fuel_mass_cruise_breguet * V_cruise / R[i]
    # print(fuel_mass_cruise_breguet)
    if pick_breguet:
        fuel_consumption.loc['cruise_breguet'] = [fuel_flow_cruise_breguet, fuel_mass_cruise_breguet]

    """
    cruise
    """
    cruise_thrust = get_thrust_required(isa(H_m)[2], V_cruise, S, CDcruise[i])/2
    # print(cruise_thrust/thrust_max)
    fuel_flow_cruise, fuel_mass_cruise = get_fuel_consumption(cruise_thrust, R[i], V_cruise)
    if not pick_breguet:
        fuel_consumption.loc['cruise'] = [fuel_flow_cruise*2, fuel_mass_cruise*2]

    if pick_breguet:
        mass -= fuel_consumption.loc['cruise_breguet']['fuel_mass']
    else:
        mass -= fuel_consumption.loc['cruise']['fuel_mass']

    # engine_failure=True
    # engines_operative = 2 - engine_failure
    """
    descent
    """
    fuel_fraction_descent = 0.990
    fuel_mass_descent = (1 - fuel_fraction_descent) * MTOW[i]
    fuel_consumption.loc['descent'] = [np.NaN, fuel_mass_descent]

    mass -= fuel_mass_descent

    """
    loiter
    """
    loiter_velocity = 130.
    loiter_altitude = 450.
    loiter_time = 5*60

    CDloiter = CDcruise[i]

    loiter_thrust = get_thrust_required(isa(loiter_altitude)[2], loiter_velocity, S, CDloiter)
    fuel_flow_loiter = get_fuel_consumption(loiter_thrust, 1, loiter_velocity)[0]
    fuel_mass_loiter = loiter_time*fuel_flow_loiter
    if pick_breguet:
        fuel_mass_trip = fuel_mass_take_off + fuel_mass_climb + fuel_mass_cruise_breguet
    else:
        fuel_mass_trip = fuel_mass_take_off + fuel_mass_climb + fuel_mass_cruise
    if fuel_mass_loiter < 0.05 * fuel_mass_trip:
        fuel_mass_loiter = 0.05 * fuel_mass_trip

    fuel_consumption.loc['loiter'] = [fuel_flow_loiter, fuel_mass_loiter]

    """
    landing
    """
    landing_field_length, landing_velocity = get_landing_field_length(engine_failure, isa(0)[2], g, h_screen_la, mass,
                                                                      thrust_max, C_L_la, C_D_la, S,
                                                                      friction_coefficient_la, reverse_thrust_factor)
    fuel_flow_landing, fuel_mass_landing = get_fuel_consumption(engines_operative * reverse_thrust_factor * thrust_max,
                                                                landing_field_length, landing_velocity/np.sqrt(2))
    fuel_consumption.loc['landing'] = [fuel_flow_landing, fuel_mass_landing]

    mass -= fuel_mass_landing

    """
    take-off #2
    """
    take_off_2_field_length, take_off_2_velocity = get_take_off_field_length(engine_failure, isa(0)[2], g, h_screen_to, mass, thrust_max, thrust_transition_setting, thrust_climb_out_setting, C_L_to, C_D_to, S, friction_coefficient_to, True, 80, reverse_thrust_factor)
    fuel_flow_take_off_2, fuel_mass_take_off_2 = get_fuel_consumption(engines_operative * thrust_max, take_off_2_field_length, take_off_2_velocity/np.sqrt(2))
    fuel_consumption.loc['take_off_2'] = [fuel_flow_take_off_2, fuel_mass_take_off_2]

    """
    climb #2
    """
    fuel_fraction_climb_2 = 0.980
    fuel_mass_climb_2 = (1 - fuel_fraction_climb_2) * MTOW[i]
    fuel_consumption.loc['climb_2'] = ['variable', fuel_mass_climb_2]
    mass -= fuel_mass_climb_2

    """
    cruise breguet #2
    """
    alternate_range = 200E3

    fuel_mass_cruise_breguet_2 = get_fuel_burned_breguet(mass, alternate_range, V_cruise, cj, g, LoverD[i])  # todo; update cj
    fuel_flow_cruise_breguet_2 = fuel_mass_cruise_breguet_2 * V_cruise / alternate_range
    if pick_breguet:
        fuel_consumption.loc['cruise_breguet_2'] = [fuel_flow_cruise_breguet_2, fuel_mass_cruise_breguet_2]

    """
    cruise #2
    """
    cruise_2_thrust = get_thrust_required(isa(H_m)[2], V_cruise, S, CDcruise[i]) / 2
    fuel_flow_cruise_2, fuel_mass_cruise_2 = get_fuel_consumption(cruise_2_thrust, alternate_range, V_cruise)
    if not pick_breguet:
        fuel_consumption.loc['cruise_2'] = [fuel_flow_cruise_2 * 2, fuel_mass_cruise_2 * 2]

    if pick_breguet:
        mass -= fuel_consumption.loc['cruise_breguet_2']['fuel_mass']
    else:
        mass -= fuel_consumption.loc['cruise_2']['fuel_mass']

    """
    descent #2
    """
    fuel_fraction_descent_2 = 0.990
    fuel_mass_descent_2 = (1 - fuel_fraction_descent_2) * MTOW[i]
    fuel_consumption.loc['descent_2'] = [np.NaN, fuel_mass_descent_2]

    mass -= fuel_mass_descent_2

    """
    loiter #2
    """
    loiter_2_time = 30 * 60
    fuel_flow_loiter_2 = get_fuel_consumption(loiter_thrust, 1, loiter_velocity)[0]
    fuel_mass_loiter_2 = loiter_2_time * fuel_flow_loiter_2
    fuel_consumption.loc['loiter_2'] = [fuel_flow_loiter_2, fuel_mass_loiter_2]

    """
    landing #2
    """
    landing_2_field_length, landing_2_velocity = get_landing_field_length(engine_failure, isa(0)[2], g, h_screen_la, mass,
                                                                      thrust_max, C_L_la, C_D_la, S,
                                                                      friction_coefficient_la, reverse_thrust_factor)
    fuel_flow_landing_2, fuel_mass_landing_2 = get_fuel_consumption(engines_operative * reverse_thrust_factor * thrust_max,
                                                                landing_2_field_length, landing_2_velocity / np.sqrt(2))
    fuel_consumption.loc['landing_2'] = [fuel_flow_landing_2, fuel_mass_landing_2]

    mass -= fuel_mass_landing_2

    return fuel_consumption


fuel_consumption = analyze_fuel_consumption(MTOW, 70)
print(sum(fuel_consumption['fuel_mass']))
fuel_consumption['fuel_mass'] = fuel_consumption['fuel_mass'] * 1.1
print(sum(fuel_consumption['fuel_mass']))

# print(fuel_consumption.loc['take_off'])

# x, y, energy_height = get_energy_height()
# plt.figure('energy height')
# plt.contour(x, y, energy_height, 30)
