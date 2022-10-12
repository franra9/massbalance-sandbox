#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:46:30 2022

@author: francesc
"""
import xarray as xr
import params as params
from params import * 
from get_data_over_glacier import get_raster_data
from projection_functions import omnibus_future
from params import *
import pandas as pd
import params as params
from totalMB import (omnibus_minimize_mf)
from scipy.optimize import minimize_scalar
import numpy as np
from get_data_over_glacier import get_raster_data
import matplotlib.pyplot as plt
import time
import warnings

# get data from rasters
n = params.n

in_data = get_raster_data(n)

# open calibration m_f raster:
o = xr.open_dataset(f'{out_path}/{y_alfa}-{y_omega}_{n}x{n}_{cal}_{ssp}_melt_f_0.nc')
in_data[2] = o.to_array()

#points onto glacier
g_ind = np.where(in_data[4].values.flatten() > 0) 

now=time.time()
warnings.filterwarnings("ignore")

# loop over each point over the glacier
melt_y = []
for i in g_ind[0]:
    
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
            in_data[1].values.flatten()[i] # geometric altitude in y_omega, on glacier sfc
    
    melt_y_i = omnibus_future(melt_f = in_data[2].values.flatten()[i], altitude = altitude, \
                   thick20 = in_data[3].values.flatten()[i], years_fcst = years_fcst)
    melt_y.append(melt_y_i)
    print(f'Point altitude is: {altitude} and its melt year is: {melt_y_i}')

now1=time.time()
print(f'Time running future evolution is: {abs(now - now1)}')

# grid size
dum = np.zeros(n * n)
dum[g_ind] = melt_y
dum_resh = np.reshape(dum, [n, n])

# fill melt_f in in_data list
in_data[2].values[0] = dum_resh

# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(13, 3), ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
pos = ax1.imshow(in_data[2].values[0], vmin=min(in_data[2].values[0].flatten()[in_data[2].values[0].flatten()>0]-10), vmax=max(in_data[2].values[0].flatten()))

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near
ax1.set_title(f'{y_ini}-{y_fin}_{n}x{n}_{cal}_{ssp}')
fig.colorbar(pos, ax=ax1)
fig.show()  
plt.savefig(f'{out_path}/{y_ini}-{y_fin}_{n}x{n}_{cal}_{ssp}_melt_year_0')
plt.show()

aaa = in_data[2]

aaa.to_netcdf(f'{out_path}/{y_ini}-{y_fin}_{n}x{n}_{cal}_{ssp}_melt_year_0.nc')
