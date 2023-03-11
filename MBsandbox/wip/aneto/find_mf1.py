#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calling functions from totalMB. To get a minimized MF using altitude change. 
(omnibus_minimize_mf1 function instead of omnibus_minimize_mf)
"""
from params import *
import params as params
from totalMB import (omnibus_minimize_mf1)
from scipy.optimize import minimize_scalar
import numpy as np
from get_data_over_glacier import get_raster_data
import matplotlib.pyplot as plt
import time
import warnings
from climate import get_climate

# get data from rasters
n = params.n

in_data = get_raster_data(n)
g_ind = np.where(in_data[4].values.flatten() > 0) #points onto glacier

now=time.time()
warnings.filterwarnings("ignore")

# loop over each point over the glacier
m_f = []
for i in g_ind[0]:
    
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
            in_data[1].values.flatten()[i]
    obs_mb = in_data[1].values.flatten()[i]
    
    obs_mb = obs_mb * 1000 # (in mm)
    obs_mb = obs_mb * rho # in mmwe

    cl = get_climate(years, altitude)    

    res = minimize_scalar(omnibus_minimize_mf1, args=(cl, altitude, obs_mb, years), tol=0.01) # check mm and mmwe in side omnibus function
    m_f.append(res.x)
    print(f'Point altitude is: {altitude} and its melt factor: {res.x}')

now1=time.time()
print(f'Time minimizing is: {abs(now - now1)}')

# grid size
dum = np.zeros(n * n)
dum[g_ind] = m_f
dum_resh = np.reshape(dum, [n, n])

# fill melt_f in in_data list
in_data[2].values = dum_resh

# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(13, 3), ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
pos = ax1.imshow(in_data[2], vmin=min(in_data[2].values.flatten()[in_data[2].values.flatten()>0]-10), vmax=max(in_data[2].values.flatten()))

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near
ax1.set_title(f'{y_alfa}-{y_omega}_{n}x{n}_{cal}_mf1_{ssp}')
fig.colorbar(pos, ax=ax1)
fig.show()  
plt.savefig(f'{out_path}/{y_alfa}-{y_omega}_{n}x{n}_{cal}_mf1_{ssp}_melt_f_0')
plt.show()

aaa = in_data[2]

aaa.to_netcdf(f'{out_path}/{y_alfa}-{y_omega}_{n}x{n}_{cal}_mf1_{ssp}_melt_f_0.nc')

