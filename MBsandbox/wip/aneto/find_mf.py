#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calling functions from totalMB. To get a minimized MF
"""
from params import *
import params as params
from scipy.optimize import minimize_scalar
import numpy as np
from get_data_over_glacier import get_raster_data
import matplotlib.pyplot as plt
import time
import warnings
import multiprocessing as mp
import sys
import os

from totalMB import omnibus_minimize_mf # extremely ugly, i know

# get data from rasters

n = params.n

in_data = get_raster_data(n)
g_ind = np.where(in_data[4].values.flatten() > 0) #points onto glacier

now=time.time()
warnings.filterwarnings("ignore")

# loop over each point over glacier
def paral(i):
    m_f = []
    #for i in g_ind[0]:
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
                in_data[1].values.flatten()[i]
    obs_mb = in_data[1].values.flatten()[i]
        
    obs_mb = obs_mb * 1000 # (in mm)
    obs_mb = obs_mb * rho # in mmwe
    
    #cl = get_climate(years, altitude)    
    
    # spinup:
    #pd_bucket0 = spinup(years, altitude)
        
    res = minimize_scalar(omnibus_minimize_mf,np.array([40,55]) ,args=(altitude, obs_mb, years), tol=0.001) # check mm and mmwe in side omnibus function
    #m_f.append(res.x) np.array([40,50])
    #    print(f'Point altitude is: {altitude} and its melt factor: {res.x}')
    return res.x #(m_f)


pool = mp.Pool(ncores) #4 cores

now=time.time()
m_f = pool.map(paral, [i for i in g_ind[0]])
now1=time.time()
print(f'Time minimizing is: {abs(now - now1)}')

pool.close()

#m_f = []
#for i in g_ind[0]:
#    
#    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
#    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
#            in_data[1].values.flatten()[i]
#    obs_mb = in_data[1].values.flatten()[i]
#    
#    obs_mb = obs_mb * 1000 # (in mm)
#    obs_mb = obs_mb * rho # in mmwe
#
#    cl = get_climate(years, altitude)    
#
#    res = minimize_scalar(omnibus_minimize_mf, args=(cl, altitude, obs_mb, years), tol=0.01) # check mm and mmwe in side omnibus function
#    m_f.append(res.x)
#    print(f'Point altitude is: {altitude} and its melt factor: {res.x}')
#
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
pos = ax1.imshow(in_data[2], vmin=min(in_data[2].values.flatten()[in_data[2].values.flatten()>0]-1), vmax=max(in_data[2].values.flatten()))

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

dir = f'{out_path}/{n}x{n}/{cal}_{ssp}/'
if not os.path.exists(dir):
    os.makedirs(dir)

ax1.set_title(f'{y_alfa}-{y_omega}+2_{n}x{n}_{cal}_{ssp}_{wspinup}{ensamble_name}_{int(rho*1000)}_none_new')
fig.colorbar(pos, ax=ax1)
#fig.show()  
plt.savefig(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}+1_{int(rho*1000)}_none_new')
#plt.show()

aaa = in_data[2]

aaa.to_netcdf(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}+1_{int(rho*1000)}_none_new.nc')

