import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Structure.Wing.isa import isa
from modules.performance.class2_performance_defs import get_take_off_field_length, get_landing_field_length,\
    get_fuel_consumption, get_climb_optimization, get_fuel_burned_breguet, get_thrust_required, get_descent
from modules.performance.serviceable_airports import *


class Performance:
    def __init__(self, C_L_to, C_L_la, C_L_cruise, C_D_0, C_D_to, C_D_la, C_D_cruise, S, OEW, MTOW, g, screen_height_to,
                 screen_height_la, thrust_max, friction_coefficient_to, friction_coefficient_la, reverse_thrust_factor,
                 engine_failure, thrust_setting_climb_out, thrust_setting_transition, payload_mass, fuel_mass_old,
                 max_airport_altitude, altitude_resolution, mass_resolution, thrust_setting_climb, altitude_cruise,
                 cruise_velocity, flying_range, lift_over_drag, aspect_ratio, oswald_efficiency_number,
                 correction_factor_to, show_plots, show_airport_plots, thrust_setting_descent):
        self.C_L_to = C_L_to
        self.C_D_to = C_D_to
        self.C_L_la = C_L_la
        self.C_D_la = C_D_la
        self.C_L_cruise = C_L_cruise
        self.C_D_cruise = C_D_cruise
        self.S = S
        self.OEW = OEW
        self.MTOW = MTOW
        self.g = g
        self.screen_height_to = screen_height_to
        self.screen_height_la = screen_height_la
        self.thrust_max = thrust_max
        self.friction_coefficient_to = friction_coefficient_to
        self.friction_coefficient_la = friction_coefficient_la
        self.reverse_thrust_factor = reverse_thrust_factor
        self.engine_failure = engine_failure
        self.thrust_setting_climb_out = thrust_setting_climb_out
        self.thrust_setting_transition = thrust_setting_transition
        self.thrust_setting_climb = thrust_setting_climb
        self.thrust_setting_descent = thrust_setting_descent
        self.payload_mass = payload_mass
        self.fuel_mass_old = fuel_mass_old
        self.max_airport_altitude = max_airport_altitude
        self.altitude_resolution = altitude_resolution
        self.mass_resolution = mass_resolution
        self.altitude_cruise = altitude_cruise
        self.cruise_velocity = cruise_velocity
        self.flying_range = flying_range
        self.lift_over_drag = lift_over_drag
        self.aspect_ratio = aspect_ratio
        self.C_D_0 = C_D_0
        self.oswald_efficiency_number = oswald_efficiency_number
        self.correction_factor_to = correction_factor_to
        self.show_plots = show_plots
        self.show_airport_plots = show_airport_plots

        self.cj = self.cj()
        self.tofl, self.airport_altitude_list, self.take_off_field_length, self.take_off_velocity, self.decision_speed = self.analyze_take_off_performance()
        self.lfl, self.landing_field_length, self.approach_velocity = self.analyze_landing_performance()
        self.fuel_table, self.fuel_mass_engine_startup, self.fuel_mass_taxi, self.fuel_mass_climb,\
        self.fuel_mass_cruise_breguet, self.fuel_mass_descent, self.fuel_mass_loiter, self.fuel_mass_landing,\
        self.fuel_mass_take_off_2, self.fuel_mass_climb_2, self.fuel_mass_cruise_breguet_2, self.fuel_mass_descent_2,\
        self.fuel_mass_loiter_2, self.fuel_mass_landing_2, self.fuel_flow_take_off, self.fuel_flow_climb,\
        self.fuel_flow_cruise_breguet, self.fuel_flow_loiter, self.fuel_flow_landing, self.fuel_flow_take_off_2,\
        self.fuel_flow_climb_2, self.fuel_flow_cruise_breguet_2, self.fuel_flow_loiter_2, self.fuel_flow_landing_2,\
        self.fuel_mass_total, self.fuel_mass_nominal, self.fuel_fraction_total, self.fuel_flow_descent,\
        self.fuel_flow_descent_2, self.fuel_mass_take_off, self.fuel_fraction_take_off, self.fuel_fraction_climb,\
        self.fuel_fraction_cruise_breguet, self.fuel_fraction_descent, self.fuel_fraction_loiter,\
        self.fuel_fraction_landing, self.fuel_fraction_take_off, self.fuel_fraction_climb_2,\
        self.fuel_fraction_cruise_breguet_2, self.fuel_fraction_descent_2, self.fuel_fraction_loiter_2,\
        self.fuel_fraction_landing_2, self.fuel_fraction_take_off_2, fuel_fraction_descent_2 = self.analyze_fuel_consumption()
        print(self.fuel_table)

    def cj(self):
        cj = get_fuel_consumption(self.thrust_max, 1, 1)[0] / self.thrust_max
        return cj

    def analyze_take_off_performance(self):
        # assumes maximum thrust

        airport_altitude_list = np.linspace(0, self.max_airport_altitude, self.altitude_resolution)
        tofl = []
        mass_list = np.linspace(self.OEW, self.MTOW, self.mass_resolution)
        # todo; review maximum/minimum take-off weight

        for altitude in airport_altitude_list:
            density = isa(altitude)[2]
            take_off_field_length_list = []; take_off_velocity_list = []
            decision_speed_list = []
            for mass in mass_list:
                length, take_off_velocity, decision_speed = get_take_off_field_length(self.engine_failure, density,
                                                                                      self.g, self.screen_height_to,
                                                                                      mass, self.thrust_max,
                                                                                      self.thrust_setting_transition,
                                                                                      self.thrust_setting_climb_out,
                                                                                      self.C_L_to, self.C_D_to, self.S,
                                                                                      self.friction_coefficient_to,
                                                                                      True, 30,
                                                                                      self.reverse_thrust_factor)
                take_off_field_length_list.append(length * self.correction_factor_to)
                take_off_velocity_list.append(take_off_velocity)
                decision_speed_list.append(decision_speed)
            tofl.append([mass_list, take_off_field_length_list, take_off_velocity_list, decision_speed_list])

        'selecting sea-level, MTOW values'
        tofl_select = tofl[0]
        take_off_field_length, take_off_velocity, decision_speed = tofl_select[1][-1], tofl_select[2][-1], tofl_select[3][-1]

        # plotting
        plt.figure('take off field length')
        for select in tofl:
            h = airport_altitude_list[tofl.index(select)]
            plt.plot(select[1], select[0], label='%a [m]' % h)

        plt.legend()
        engines_used = 2 - self.engine_failure

        plt.title('Take-off field length, %a engine(s) operative' % engines_used)
        plt.axhline(y=self.OEW, linestyle=':')
        plt.text(2000, self.OEW, 'OEW')
        plt.axhline(y=self.MTOW, linestyle=':')
        plt.text(2000, self.MTOW, 'MTOW')
        plt.axhline(y=self.OEW + self.payload_mass, linestyle=':')
        plt.text(2000, self.OEW + self.payload_mass, 'OEW+M_payload')
        plt.axhline(y=self.OEW + self.fuel_mass_old, linestyle=':')
        plt.text(2000, self.OEW + self.fuel_mass_old, 'OEW+M_fuel')

        plt.axvline(2000, linestyle=':')

        plt.ylabel('Mass [kg]')
        plt.xlabel('Take-off field length [m]')
        if self.show_plots is False:
            plt.close()

        return tofl, airport_altitude_list, take_off_field_length, take_off_velocity, decision_speed

    def analyze_landing_performance(self):
        correction_factor_la = 1.   # todo; review

        airport_altitude_list = np.linspace(0, self.max_airport_altitude, self.altitude_resolution)
        lfl = []
        mass_list = np.linspace(self.OEW, self.MTOW, self.mass_resolution)
        # todo; review minimum/maximum landing weight

        for altitude in airport_altitude_list:
            density = isa(altitude)[2]
            landing_field_length_list = []
            approach_velocity_list = []
            for mass in mass_list:
                x, approach_velocity = get_landing_field_length(self.engine_failure, density, self.g,
                                                                self.screen_height_la, mass, self.thrust_max,
                                                                self.C_L_la, self.C_D_la, self.S,
                                                                self.friction_coefficient_la,
                                                                self.reverse_thrust_factor)
                landing_field_length_list.append(x * correction_factor_la)
                approach_velocity_list.append(approach_velocity)
            lfl.append([mass_list, landing_field_length_list, approach_velocity_list])

        'selecting sea-level, MTOW values'
        lfl_select = lfl[0]
        landing_field_length, approach_velocity = lfl_select[1][-1], lfl_select[2][-1]

        # plotting
        plt.figure('landing field length')
        for select in lfl:
            h = airport_altitude_list[lfl.index(select)]
            plt.plot(select[1], select[0], label='%a [m]' % h)

        plt.legend()
        engines_used = 2 - self.engine_failure
        plt.axvline(1300, linestyle=':')
        plt.title('Landing field length, %a engine(s) operative' % engines_used)
        # plt.plot([0, max(select[1])], [MTOW[i], MTOW[i]], color='C0', linestyle=':')
        # plt.plot([0, max(select[1])], [OEW[i], OEW[i]], color='C0', linestyle=':')
        plt.ylabel('Mass [kg]')
        plt.xlabel('Landing field length [m]')
        if self.show_plots is False:
            plt.close()
        return lfl, landing_field_length, approach_velocity

    def analyze_fuel_consumption(self):
        """
        inputs
        """
        pick_breguet = True
        alternate_range = 200E3         # todo; review
        loiter_velocity = 130.          # todo; review
        loiter_altitude = 450.          # ICAO [s]
        loiter_time = 5 * 60            # ICAO [s]
        loiter_2_time = 30 * 60         # ICAO [s]
        altitude_cruise_2 = self.altitude_cruise * 0.5

        """
        analysis
        """
        fuel_consumption = pd.DataFrame(data=None, columns=['fuel_flow', 'fuel_mass'])

        engines_operative = 2 - self.engine_failure

        'engine_startup, taxi'
        fuel_fraction_engine_startup = 0.990
        fuel_fraction_taxi = 0.990

        fuel_mass_engine_startup = (1 - fuel_fraction_engine_startup) * self.MTOW
        fuel_mass_taxi = (1 - fuel_fraction_taxi) * self.MTOW
        fuel_consumption.loc['engine_startup'] = [np.NaN, fuel_mass_engine_startup]
        fuel_consumption.loc['taxi'] = [np.NaN, fuel_mass_taxi]

        'take_off'
        # sea level
        # assume full thrust during entire take-off
        take_off_field_length, take_off_velocity = get_take_off_field_length(self.engine_failure, isa(0)[2], self.g,
                                                                             self.screen_height_to, self.MTOW,
                                                                             self.thrust_max,
                                                                             self.thrust_setting_transition,
                                                                             self.thrust_setting_climb_out, self.C_L_to,
                                                                             self.C_D_to, self.S,
                                                                             self.friction_coefficient_to, True, 30,
                                                                             self.reverse_thrust_factor)[:2]
        fuel_flow_take_off, fuel_mass_take_off = get_fuel_consumption(engines_operative * self.thrust_max,
                                                                      take_off_field_length,
                                                                      take_off_velocity / np.sqrt(2))
        fuel_consumption.loc['take_off'] = [fuel_flow_take_off, fuel_mass_take_off]

        mass = self.MTOW - fuel_consumption.loc['take_off']['fuel_mass'] - fuel_consumption.loc['engine_startup'][
            'fuel_mass'] - fuel_consumption.loc['taxi']['fuel_mass']

        'climb'
        fuel_flow_climb, fuel_mass_climb, climb_final_velocity, distance_climb = get_climb_optimization(mass, self.thrust_max,
                                                                                                  self.C_D_cruise,
                                                                                                  self.S, self.g,
                                                                                                  self.altitude_cruise,
                                                                                                  self.cruise_velocity,
                                                                                                  self.thrust_setting_climb)
        fuel_consumption.loc['climb'] = [fuel_flow_climb, fuel_mass_climb]
        mass -= fuel_mass_climb
        if self.show_plots is False:
            plt.close()
        'cruise_breguet'
        range_leftover = self.flying_range - distance_climb
        fuel_mass_cruise_breguet = get_fuel_burned_breguet(mass, range_leftover, self.cruise_velocity, self.cj, self.g,
                                                           self.lift_over_drag)
        fuel_flow_cruise_breguet = fuel_mass_cruise_breguet * self.cruise_velocity / range_leftover

        if pick_breguet:
            fuel_consumption.loc['cruise_breguet'] = [fuel_flow_cruise_breguet, fuel_mass_cruise_breguet]

        'cruise_non_breguet'
        cruise_thrust = get_thrust_required(isa(self.altitude_cruise)[2], self.cruise_velocity, self.S, self.C_D_cruise) / 2
        fuel_flow_cruise, fuel_mass_cruise = get_fuel_consumption(cruise_thrust, range_leftover, self.cruise_velocity)

        if pick_breguet:
            mass -= fuel_consumption.loc['cruise_breguet']['fuel_mass']
        else:
            mass -= fuel_consumption.loc['cruise']['fuel_mass']
            fuel_consumption.loc['cruise'] = [fuel_flow_cruise * 2, fuel_mass_cruise * 2]

        'descent'
        fuel_flow_descent = get_fuel_consumption(self.thrust_setting_descent * self.thrust_max, 1, 1)[0]
        fuel_mass_descent, distance_descent = get_descent(self.altitude_cruise, self.cruise_velocity,
                                                          self.approach_velocity, engines_operative,
                                                          self.thrust_setting_descent * self.thrust_max,
                                                          self.S, self.C_D_cruise, mass, self.g)
        fuel_consumption.loc['descent'] = [fuel_flow_descent, fuel_mass_descent]
        mass -= fuel_mass_descent

        'loiter'
        loiter_thrust = get_thrust_required(isa(loiter_altitude)[2], loiter_velocity, self.S, self.C_D_cruise)
        fuel_flow_loiter = get_fuel_consumption(loiter_thrust, 1, loiter_velocity)[0]
        fuel_mass_loiter = loiter_time * fuel_flow_loiter
        if pick_breguet:
            fuel_mass_trip = fuel_mass_take_off + fuel_mass_climb + fuel_mass_cruise_breguet
        else:
            fuel_mass_trip = fuel_mass_take_off + fuel_mass_climb + fuel_mass_cruise
        if fuel_mass_loiter < 0.05 * fuel_mass_trip:
            fuel_mass_loiter = 0.05 * fuel_mass_trip

        fuel_consumption.loc['loiter'] = [fuel_flow_loiter, fuel_mass_loiter]

        'landing'
        landing_field_length, landing_velocity = get_landing_field_length(self.engine_failure, isa(0)[2], self.g,
                                                                          self.screen_height_la, mass, self.thrust_max,
                                                                          self.C_L_la, self.C_D_la, self.S,
                                                                          self.friction_coefficient_la,
                                                                          self.reverse_thrust_factor)
        fuel_flow_landing, fuel_mass_landing = get_fuel_consumption(
            engines_operative * self.reverse_thrust_factor * self.thrust_max,
            landing_field_length, landing_velocity / np.sqrt(2))
        fuel_consumption.loc['landing'] = [fuel_flow_landing, fuel_mass_landing]

        mass -= fuel_mass_landing

        'take_off_2'
        take_off_2_field_length, take_off_2_velocity = get_take_off_field_length(self.engine_failure, isa(0)[2], self.g,
                                                                                 self.screen_height_to, mass,
                                                                                 self.thrust_max,
                                                                                 self.thrust_setting_transition,
                                                                                 self.thrust_setting_climb_out,
                                                                                 self.C_L_to, self.C_D_to, self.S,
                                                                                 self.friction_coefficient_to,
                                                                                 True, 30, self.reverse_thrust_factor)[:2]
        fuel_flow_take_off_2, fuel_mass_take_off_2 = get_fuel_consumption(engines_operative * self.thrust_max,
                                                                          take_off_2_field_length,
                                                                          take_off_2_velocity / np.sqrt(2))
        fuel_consumption.loc['take_off_2'] = [fuel_flow_take_off_2, fuel_mass_take_off_2]

        'climb_2'
        fuel_flow_climb_2, fuel_mass_climb_2, climb_2_final_velocity, distance_climb_2 = \
            get_climb_optimization(mass, self.thrust_max, self.C_D_cruise, self.S, self.g, altitude_cruise_2,
                                   self.cruise_velocity, self.thrust_setting_climb)
        fuel_consumption.loc['climb_2'] = [fuel_flow_climb_2, fuel_mass_climb_2]
        mass -= fuel_mass_climb_2
        if self.show_plots is False:
            plt.close()

        'cruise_breguet_2'
        fuel_mass_cruise_breguet_2 = get_fuel_burned_breguet(mass, alternate_range, self.cruise_velocity, self.cj,
                                                             self.g, self.lift_over_drag)
        fuel_flow_cruise_breguet_2 = fuel_mass_cruise_breguet_2 * self.cruise_velocity / alternate_range

        if pick_breguet:
            fuel_consumption.loc['cruise_breguet_2'] = [fuel_flow_cruise_breguet_2, fuel_mass_cruise_breguet_2]

        'cruise_non_breguet_2'
        cruise_2_thrust = get_thrust_required(isa(altitude_cruise_2)[2], self.cruise_velocity, self.S, self.C_D_cruise) / 2
        fuel_flow_cruise_2, fuel_mass_cruise_2 = get_fuel_consumption(cruise_2_thrust, alternate_range,
                                                                      self.cruise_velocity)

        if pick_breguet:
            mass -= fuel_consumption.loc['cruise_breguet_2']['fuel_mass']
        else:
            mass -= fuel_consumption.loc['cruise_2']['fuel_mass']
            fuel_consumption.loc['cruise_2'] = [fuel_flow_cruise_2 * 2, fuel_mass_cruise_2 * 2]

        'descent_2'
        fuel_flow_descent_2 = get_fuel_consumption(self.thrust_setting_descent * self.thrust_max, 1, 1)[0]
        fuel_mass_descent_2, distance_descent_2 = get_descent(altitude_cruise_2, self.cruise_velocity,
                                                          self.approach_velocity, engines_operative,
                                                          self.thrust_setting_descent * self.thrust_max,
                                                          self.S, self.C_D_cruise, mass, self.g)
        fuel_consumption.loc['descent_2'] = [fuel_flow_descent_2, fuel_mass_descent_2]

        mass -= fuel_mass_descent_2

        'loiter_2'
        fuel_flow_loiter_2 = get_fuel_consumption(loiter_thrust, 1, loiter_velocity)[0]
        fuel_mass_loiter_2 = loiter_2_time * fuel_flow_loiter_2
        fuel_consumption.loc['loiter_2'] = [fuel_flow_loiter_2, fuel_mass_loiter_2]

        'landing_2'
        landing_2_field_length, landing_2_velocity = get_landing_field_length(self.engine_failure, isa(0)[2], self.g,
                                                                              self.screen_height_la, mass,
                                                                              self.thrust_max, self.C_L_la, self.C_D_la,
                                                                              self.S, self.friction_coefficient_la,
                                                                              self.reverse_thrust_factor)
        fuel_flow_landing_2, fuel_mass_landing_2 = get_fuel_consumption(
            engines_operative * self.reverse_thrust_factor * self.thrust_max,
            landing_2_field_length, landing_2_velocity / np.sqrt(2))
        fuel_consumption.loc['landing_2'] = [fuel_flow_landing_2, fuel_mass_landing_2]

        mass -= fuel_mass_landing_2
        fuel_mass_total = sum(fuel_consumption['fuel_mass'])
        fuel_mass_nominal = sum([fuel_mass_engine_startup, fuel_mass_take_off, fuel_mass_climb, fuel_mass_cruise_breguet, fuel_mass_descent, fuel_mass_landing])
        fuel_fraction_total = 1 - fuel_mass_total/self.MTOW
        fuel_consumption['fuel_fraction'] = 1 - fuel_consumption['fuel_mass'] / self.MTOW
        fuel_fraction_take_off = 1 - fuel_mass_take_off / self.MTOW
        fuel_fraction_climb = 1 - fuel_mass_climb / self.MTOW
        fuel_fraction_cruise_breguet = 1 - fuel_mass_cruise_breguet / self.MTOW
        fuel_fraction_loiter = 1 - fuel_mass_loiter / self.MTOW
        fuel_fraction_descent = 1 - fuel_mass_descent / self.MTOW
        fuel_fraction_landing = 1 - fuel_mass_landing / self.MTOW
        fuel_fraction_take_off_2 = 1 - fuel_mass_take_off_2 / self.MTOW
        fuel_fraction_climb_2 = 1 - fuel_mass_climb_2 / self.MTOW
        fuel_fraction_cruise_breguet_2 = 1 - fuel_mass_cruise_breguet_2 / self.MTOW
        fuel_fraction_loiter_2 = 1 - fuel_mass_loiter_2 / self.MTOW
        fuel_fraction_landing_2 = 1 - fuel_mass_landing_2 / self.MTOW
        fuel_fraction_descent_2 = 1 - fuel_mass_descent_2 / self.MTOW
        return fuel_consumption, fuel_mass_engine_startup, fuel_mass_taxi, fuel_mass_climb, fuel_mass_cruise_breguet, fuel_mass_descent, fuel_mass_loiter, fuel_mass_landing, fuel_mass_take_off_2, fuel_mass_climb_2, fuel_mass_cruise_breguet_2, fuel_mass_descent_2, fuel_mass_loiter_2, fuel_mass_landing_2, fuel_flow_take_off, fuel_flow_climb, fuel_flow_cruise_breguet, fuel_flow_loiter, fuel_flow_landing, fuel_flow_take_off_2, fuel_flow_climb_2, fuel_flow_cruise_breguet_2, fuel_flow_loiter_2, fuel_flow_landing_2, fuel_mass_total, fuel_mass_nominal, fuel_fraction_total, fuel_flow_descent, fuel_flow_descent_2, fuel_mass_take_off, fuel_fraction_take_off, fuel_fraction_climb, fuel_fraction_cruise_breguet, fuel_fraction_descent, fuel_fraction_loiter, fuel_fraction_landing, fuel_fraction_take_off, fuel_fraction_climb_2, fuel_fraction_cruise_breguet_2, fuel_fraction_descent_2, fuel_fraction_loiter_2, fuel_fraction_landing_2, fuel_fraction_take_off_2, fuel_fraction_descent_2

    def get_serviceable_airports(self):
        serviceable_airports(self.landing_field_length, self.airport_altitude_list, self.flying_range, self.show_airport_plots)

# C_L_to = 1.9
# C_L_la = 2.3
# C_D_0 = 0.018117539865047032
# C_D_to = 0.19542078310015343
# C_D_la = 0.28017202214599246
# S = 132.
# OEW = 34631.92
# MTOW = 59387.234177062084
# g = 9.81
# thrust_max = thrust_max
# friction_coefficient_to = 0.05
# friction_coefficient_la = 0.45
# reverse_thrust_factor = 0.45
# thrust_setting_climb_out = 1.
# thrust_setting_transition = 1.
# payload_mass = M_payload[0]
# fuel_mass = M_fuel[0]
# thrust_setting_climb = 0.9
# altitude_cruise = H_m
# cruise_velocity = V_cruise
# flying_range = 4000000.0
# cj = 0.790/thrust_max
# lift_over_drag = LoverD[0]
# aspect_ratio = A
# oswald_efficiency_number = e
# C_L_cruise = 0.8
# C_D_cruise = 0.05

# 'inputs'
# max_airport_altitude = 2000.  # [m]
# altitude_resolution = 5  # number of different altitudes considered
# mass_resolution = 20  # resolution of plotting mass vs take-off field length
# engine_failure = False
# screen_height_to = 35 * ft_to_m
# screen_height_la = 50 * ft_to_m
#
#
#
# 'analysis'
# config1_Performance = Performance(C_L_to, C_L_la, C_L_cruise, C_D_0, C_D_to, C_D_la, C_D_cruise, S, OEW, MTOW, g, screen_height_to,
#                  screen_height_la, thrust_max, friction_coefficient_to, friction_coefficient_la, reverse_thrust_factor,
#                  engine_failure, thrust_setting_climb_out, thrust_setting_transition, payload_mass, fuel_mass,
#                  max_airport_altitude, altitude_resolution, mass_resolution, thrust_setting_climb, altitude_cruise,
#                  cruise_velocity, flying_range, cj, lift_over_drag, aspect_ratio, oswald_efficiency_number)
