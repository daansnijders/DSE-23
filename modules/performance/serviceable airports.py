import pandas as pd
import numpy as np
from inputs.concept_1 import ft_to_m
import matplotlib.pyplot as plt
from modules.performance.class2_performance_tests import analyze_take_off_performance, i
from inputs.concept_1 import R
from matplotlib.lines import Line2D

# source; http://ourairports.com/
# We'd love you to give us credit, like we give credit to our sources, but you're not required to.

# inputs
flying_range = R[i]
altitude_length_list, airport_altitude_list = analyze_take_off_performance()

df_airports = pd.read_csv('../airports.csv', index_col=0, header=0, names=['id', 'ident', 'type', 'name', 'latitude', 'longitude',
                                                                 'elevation', 'continent', 'iso_country', 'iso_region',
                                                                 'municipality', 'scheduled_service', 'gps_code',
                                                                 'iata_code', 'local_code', 'home-link', 'wikipedia',
                                                                 'keywords'])

# filtering, cleaning)
df_airports = df_airports[pd.notnull(df_airports['elevation'])]

df_airports['elevation'] = df_airports['elevation'] * ft_to_m

df_runways = pd.read_csv('../runways.csv', index_col=0, header=0)
df_runways = df_runways[['airport_ident', 'length_ft']]
df_runways['ident'] = df_runways['airport_ident']
df_runways = df_runways[['ident', 'length_ft']]
df_runways['length_m'] = df_runways['length_ft'] * ft_to_m


df_airports = df_airports[['ident', 'type', 'name', 'latitude', 'longitude', 'elevation', 'iso_country', 'continent']]

df = pd.merge(df_airports, df_runways, on='ident')

"""
Filtering options
---
Can only look at flights in USA, only large airports, etc
"""

# airport type selection
identifiers = ['small_airport', 'medium_airport', 'large_airport']
df = df[df['type'].isin(identifiers)]

# # country selection
# identifiers = ['NL']
# df = df[df['iso_country'].isin(identifiers)]

# continent selection
# identifiers = ['EU'] # change this to check only for large airports
# df = df[df['continent'].isin(identifiers)]


plt.close()
plt.close()
plt.close()
plt.figure()

plt.scatter(df['longitude'], df['latitude'], s=1, color='C0')

# m = folium.Map(
#     location=[46.3014, -123.7390],
#     zoom_start=7,
#     tiles='Stamen Terrain'
# )

# analysis
number_unique_airports_total = len(df['name'].drop_duplicates())

servable_airports_total = pd.DataFrame()
for select in altitude_length_list:
    runway_length_MTOW = max(select[1])
    h = airport_altitude_list[altitude_length_list.index(select)]
    plt.plot(select[1], select[0], label='%a [m]' % h)

    servable_airports = df[df['elevation'] < h]
    servable_airports = servable_airports[servable_airports['length_m'] > runway_length_MTOW]
    plt.scatter(servable_airports['longitude'], servable_airports['latitude'], marker="^", s=1, color='C1')
    # for a in range(0, len(servable_airports), 1):
    #     folium.Marker([servable_airports.iloc[a]['longitude'], servable_airports.iloc[a]['latitude']], popup='<i>Mt. Hood Meadows</i>', tooltip='Serviceable').add_to(m)
    servable_airports_total = servable_airports_total.append(servable_airports)

unique_servable_airports = servable_airports_total['name'].drop_duplicates()
number_unique_servable_airports = len(unique_servable_airports)
# print(number_airports_servable/number_airports_total)
# print(servable_airports[['name', 'length_m', 'elevation']])

# plotting range around schiphol'])
longitude = (df.loc[df['ident'] == 'EHAM']).iloc[0]['longitude']
latitude = (df.loc[df['ident'] == 'EHAM']).iloc[0]['latitude']
steps = 500

circle1 = []
circle2 = []

for a in np.linspace(-1, 1, steps):
    y = a * flying_range
    x = np.sqrt(flying_range**2-y**2)
    y_deg = latitude + y/111000
    x_deg_neg = longitude - x/ (np.cos(y_deg/180*np.pi)* 111321)
    circle1.append((x_deg_neg, y_deg))

for a in np.linspace(-.999999999, .999999999, steps):
    y = a * flying_range
    x = np.sqrt(flying_range ** 2 - y ** 2)
    y_deg = latitude + y / 111000
    x_deg_pos = longitude + x / (np.cos(y_deg / 180 * np.pi) * 111321)
    circle2.append((x_deg_pos, y_deg))

testlist1 = [(elem1, elem2) for elem1, elem2 in circle1]
testlist2 = [(elem1, elem2) for elem1, elem2 in circle2]
plt.plot(*zip(*testlist1), color='C2')
plt.plot(*zip(*testlist2), color='C2')

# circle1 = Ellipse(((df.loc[df['ident'] == 'EHAM']).iloc[0]['longitude'], (df.loc[df['ident'] == 'EHAM']).iloc[0]['latitude']), R[i]*2*,  color='C2')

plt.xlim((-180, 180))
plt.ylim((-90, 90))
fig = plt.gcf()
ax = fig.gca()

legend_elements = [Line2D([0], [0], marker='o', color='C0', label='Non-serviceable airport'),
                   Line2D([0], [0], marker='o', color='C1', label='Serviceable airport')]

# Create the figure
plt.legend(handles=legend_elements)


plt.show()
print(number_unique_airports_total)
print(number_unique_servable_airports/number_unique_airports_total)

# m.save('index.html')