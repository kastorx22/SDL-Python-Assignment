'''
need to install the following first:
pip install geodatasets
pip install geopandas

'''
#Kateryna
import numpy as np
import matplotlib.pyplot as plt
#dictionary with each value
ascend_dict={
    'time':[],
    'dist':[],
    'long':[],
    'lat':[],
    'vel_az':[],
    'vel_el':[],
    'ef_vel':[],
    'head':[],
    'flt_path':[],
    'sf_vel':[],
    'range':[],
    'altitude':[]
}
keys=list(ascend_dict.keys())
#print(keys)
with open('as-505-ascent-phase-data.txt', mode='r') as file:
    for line in file:
        line_str=line.strip()
        #to skip blank lines, comments, that start with $, and headers
        if(
            not line_str
            or line_str.startswith('$')
            or line_str.startswith('TIME')
            or line_str.startswith('SEC')
        ):
            continue
        #replace comas with spaces and then skip empty strings
        row=line_str.replace(',' , ' ').split()
        #if len(row)==len(ascend_dict): has extra 0s that mess the plots
        #can avoid extra 0s by using 12 columns exactly
        if len(row)==12:
            #need to convert to float first, to avoid the error "ValueError: invalid literal for int() with base 10: '0.0'"
            int_values=[int(float(val)) for val in row]
            #assign each value to appropriate key
            for key, val in zip(ascend_dict.keys(), int_values):
                ascend_dict[key].append(val)

#convert lists to numpy array
for key in ascend_dict:
    ascend_dict[key]=np.array(ascend_dict[key])

#print(ascend_dict)
#dictionary with units for the labels
units={
    'dist':'km',
    'long': 'deg E',
    'lat':'deg N',
    'vel_az':'deg',
    'vel_el': 'deg',
    'ef_vel': 'm/s',
    'head': 'deg',
    'flt_path':'deg',
    'sf_vel':'m/s',
    'range':'m',
    'altitude':'m'
}
#parameters excluding time
parameter=[key for key in ascend_dict.keys() if key!='time']

fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12,10))
axes = axes.ravel()

for i, (col, ax) in enumerate(zip(parameter, axes)):
    #to use a unit correscopnding to each parameter
    unit=units.get(col, '')
    ax.plot( ascend_dict['time'], ascend_dict[col], label=f"{i}, {parameter}", color='m')
    ax.set_xlabel('time (s)')
    ax.set_ylabel(f"{col} ({unit})")
    ax.set_title(f"{col} vs time")
#hide the empty 12th subplot
fig.delaxes(axes[11])

plt.tight_layout()
plt.show()

#done by  Enrico
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import geodatasets
import geopandas as gpd
path = geodatasets.get_path("naturalearth.land")
world = gpd.read_file(path)

fig, ax = plt.subplots(figsize=(12, 6))
world.plot(color="lightgrey", ax=ax)

x = ascend_dict['long']
y = ascend_dict['lat']
z = ascend_dict['altitude']

plt.scatter(x, y, c=z, cmap='plasma')


plt.colorbar(label='Altitude [meter]')
plt.title("NASA: Apollo 10 ground control")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()
