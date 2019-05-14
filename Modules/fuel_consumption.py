from ISA import isa_meters


def fuel_flow(thrust):
    # only true for PW1525G
    # approximately linear relation between thrust and fuel consumption
    # from linear regression;
    ff = 7E-6 * thrust + 0.0147
    # R^2 = 0.9981
    return ff


def cruise_thrust(cruise_altitude, V_cruise, S_wing, C_D_cruise):
    T = 0.5 * isa_meters(cruise_altitude)[2] * V_cruise**2 * S_wing * C_D_cruise
    return T


def cruise_fuel(thrust_cruise): # todo fuel per kilometer (?)
    fuel_flow_cruise = fuel_flow(thrust_cruise)
    return fuel_flow_cruise