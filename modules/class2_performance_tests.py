from Structure.Wing.isa import isa
from inputs.concept_1 import ft_to_m, OEW, MTOW, thrust_max, S, P_nw, x_mlg, x_nlg, z_mlg, x_cg, z_cg, g, M_payload, M_fuel
from modules.class2_performance_defs import *
import matplotlib.pyplot as plt

"""
inputs
"""
i = 2   # configuration selection
h_screen_to = 35 * ft_to_m                                  # [m]
h_screen_la = 50 * ft_to_m


"""
take-off performance
"""


def analyze_take_off_performance():
    # todo; determine decision speed (need brake design)

    # inputs
    engine_failure = False
    thrust_transition_setting = 1.
    thrust_climb_out_setting = 1.                           # todo; reconsider this value (sometimes 0.85)
    C_L_to = 1.9                                            # todo; lift and drag of canard (update C_L_TO, C_D_TO, S)
    C_D_to = 0.001                                          # todo; lift and drag of canard (update C_L_TO, C_D_TO, S)
    correction_factor = 1.15                                # from CS-25
    # friction_coefficient = get_friction_coefficient(P_nw[i], mass, x_mlg[i], x_nlg[i], x_cg[i], z_cg - z_mlg[i], g)
    friction_coefficient = 0.05                             # from Raymer table 17.1

    max_airport_altitude = 1000.                            # [m]
    altitude_resolution = 5
    airport_altitude_list = np.linspace(0, max_airport_altitude, altitude_resolution)

    take_off_field_length = []

    mass_resolution = 50
    mass_list = np.linspace(OEW[i], MTOW[i], mass_resolution)   # todo; incorrect

    for altitude in airport_altitude_list:
        density = isa(altitude)[2]
        take_off_field_length_list = []
        for mass in mass_list:

            take_off_field_length_list.append(get_take_off_field_length(engine_failure, density, g, h_screen_to, mass,
                                                                        thrust_max, thrust_transition_setting,
                                                                        thrust_climb_out_setting, C_L_to, C_D_to, S,
                                                                        friction_coefficient) *
                                              correction_factor)
        take_off_field_length.append([mass_list, take_off_field_length_list])

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


"""
landing performance
"""


def analyze_landing_performance():
    # inputs
    engine_failure = False
    C_L_la = 2.4                    # todo; lift and drag of canard (update C_L_LA, C_D_LA, S)
    C_D_la = 0.02                  # todo; lift and drag of canard (update C_L_LA, C_D_LA, S)
    # friction_coefficient = get_friction_coefficient(P_nw[i], mass, x_mlg[i], x_nlg[i], x_cg[i], z_cg - z_mlg[i], g)
    friction_coefficient = 0.45      # from Raymer table 17.1
    correction_factor = 1.          # todo; find in CS-25??? look it up for landing
    reverse_thrust_factor = 0.45

    max_airport_altitude = 1000.  # [m]
    altitude_resolution = 5
    airport_altitude_list = np.linspace(0, max_airport_altitude, altitude_resolution)

    landing_field_length = []

    mass_resolution = 50
    mass_list = np.linspace(OEW[i], MTOW[i], mass_resolution)   # todo; review this, it's incorrect

    for altitude in airport_altitude_list:
        density = isa(altitude)[2]
        landing_field_length_list = []
        for mass in mass_list:

            landing_field_length_list.append(get_landing_field_length(engine_failure, density, g, h_screen_la, mass,
                                                                      thrust_max,C_L_la, C_D_la, S,
                                                                      friction_coefficient, reverse_thrust_factor) *
                                             correction_factor)
        landing_field_length.append([mass_list, landing_field_length_list])

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


analyze_take_off_performance()
# analyze_landing_performance()
