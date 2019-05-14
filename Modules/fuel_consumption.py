from ISA import isa_meters


def fuel_flow(thrust):
    # only true for PW1525G
    # approximately linear relation between thrust and fuel consumption
    # from linear regression;
    ff = 7E-6 * thrust + 0.0147                                                             # kg/s
    # R^2 = 0.9981
    return ff


def cruise_thrust(cruise_altitude, cruise_velocity, S_wing, C_D_cruise):
    T = 0.5 * isa_meters(cruise_altitude)[2] * cruise_velocity**2 * S_wing * C_D_cruise     # N
    return T


def cruise_fuel(thrust_cruise, cruise_range, cruise_velocity): # todo fuel per kilometer (?)
    fuel_flow_cruise = fuel_flow(thrust_cruise)
    total_fuel_used_cruise = cruise_range / cruise_velocity * fuel_flow_cruise
    return fuel_flow_cruise
