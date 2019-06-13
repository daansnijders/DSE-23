from Structure.Wing.isa import isa
from inputs.concept_1 import ft_to_m, OEW, MTOW, thrust_max, S, P_nw, x_mlg, x_nlg, z_mlg, x_cg, z_cg, g
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
            # friction_coefficient = get_friction_coefficient(P_nw[i], mass, x_mlg[i], x_nlg[i], x_cg[i], z_cg - z_mlg[i], g)
            friction_coefficient = 0.05
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
    plt.title('Take-off field, %a engine(s) operative' % engines_used)
    # plt.plot([0, max(select[1])], [MTOW[i], MTOW[i]], color='C0', linestyle=':')
    # plt.plot([0, max(select[1])], [OEW[i], OEW[i]], color='C0', linestyle=':')
    plt.ylabel('Mass [kg]')
    plt.xlabel('Take-off field length [m]')


"""
landing performance
"""


def analyze_landing_performance():
    # inputs
    engine_failure = False

    C_L_la = 2.4  # todo; lift and drag of canard (update C_L_TO, C_D_TO, S)
    C_D_la = 0.001  # todo; lift and drag of canard (update C_L_TO, C_D_TO, S)

    # correction_factor = 1.5  # todo; find in CS-25??? look it up for landing

    max_airport_altitude = 1000.  # [m]
    altitude_resolution = 5
    airport_altitude_list = np.linspace(0, max_airport_altitude, altitude_resolution)

    landing_field_length = []

    mass_resolution = 50
    mass_list = np.linspace(OEW[i], MTOW[i], mass_resolution) # todo; review this, it's incorrect
