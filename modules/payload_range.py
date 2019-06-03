import matplotlib.pyplot as plt
from modules.initialsizing_weights import *
import numpy as np


def range_breguet(M_ff, V, cj, LD_cruise, g):
    range_cruise = V/cj/g*LD_cruise*np.log(1/M_ff)
    return range_cruise


def generate_payload_range_diagram(payload_mass, fuel_mass, MTOW, design_range, cruise_velocity, cj, LD_cruise, g, OEW, i, y_lim_top):
    """
    determining reserve fuel mass, usable fuel mass, zero payload range
    """

    w11w7 = 0.992*0.990*0.9876*0.9911
    reserve_fuel_mass = (1 - w11w7) * MTOW

    fuel_mass = fuel_mass - reserve_fuel_mass

    # generate point for zero payload
    zero_payload_fuel_fraction = 1-fuel_mass/(MTOW-payload_mass)
    zero_payload_range = range_breguet(zero_payload_fuel_fraction, cruise_velocity, cj, LD_cruise, g)
    # currently missing; reserve fuel part

    # fuel_fraction = 1 - fuel_mass / (MTOW)
    # design_range = range_breguet(fuel_fraction, cruise_velocity, cj, LD_cruise, g)

    """
    plotting
    """
    plt.figure()

    x_oew = [0, design_range / 1000, zero_payload_range / 1000]
    y_oew = [OEW, OEW, OEW]

    x_reserve_fuel = [0, design_range / 1000, zero_payload_range / 1000]
    y_reserve_fuel = [OEW + reserve_fuel_mass, OEW + reserve_fuel_mass, OEW + reserve_fuel_mass]

    x_payload = [0, design_range / 1000, zero_payload_range / 1000]
    y_payload = [OEW + reserve_fuel_mass + payload_mass, OEW + reserve_fuel_mass + payload_mass, OEW + reserve_fuel_mass]

    x_total_mass = [0, design_range/1000, zero_payload_range/1000]
    y_total_mass = [OEW + reserve_fuel_mass + payload_mass, MTOW, MTOW - payload_mass]

    plt.plot(x_oew, y_oew, 'C7', label='OEW')
    plt.plot(x_oew, y_oew, 'k')
    plt.plot(x_reserve_fuel, y_reserve_fuel, 'C8', label='Reserve fuel')
    plt.plot(x_reserve_fuel, y_reserve_fuel, 'k')
    plt.plot(x_payload, y_payload, 'C0', label='Payload')
    plt.plot(x_payload, y_payload, 'k')
    plt.plot(x_total_mass, y_total_mass, 'C1', label='Fuel')
    plt.plot(x_total_mass, y_total_mass, 'k')
    plt.plot([zero_payload_range/1000, zero_payload_range/1000], [MTOW - payload_mass, 0], 'k')  # vertical end line

    # coloring in
    plt.fill_between(x_oew, y_oew, 0, color='C7')
    plt.fill_between(x_payload, y_reserve_fuel, y_oew, color='C8')
    plt.fill_between(x_payload, y_payload, y_reserve_fuel, color='C0')
    plt.fill_between(x_payload, y_total_mass, y_payload, color='C1')

    # legend
    plt.legend(fontsize=15)

    plt.xlim((0))
    plt.ylim((0, y_lim_top))
    plt.xlabel('Range [km]', fontsize=15)
    plt.ylabel('Mass [kg]', fontsize=15)
    num = i+1
    plt.title('Payload-Range Diagram Config. %a' % num, fontsize=20)

    plt.show()
