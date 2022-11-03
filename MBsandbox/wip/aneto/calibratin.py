#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 2.

Perform calibration checks using w5e5 and ensemble data.

@author: francesc
"""
import xarray as xr
import params as params
from get_data_over_glacier import get_raster_data
from projection_functions import omnibus_future
from params import *
import pandas as pd
from totalMB import (omnibus_minimize_mf)
from scipy.optimize import minimize_scalar
import numpy as np
from get_data_over_glacier import get_raster_data
import matplotlib.pyplot as plt
import time
import warnings
from totalMB import spinup, monthly_mb_sd
from climate import get_climate, get_climate_cal
import multiprocessing as mp


# get data from rasters
n = params.n

in_data = get_raster_data(n)

# open calibration m_f raster:
o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2018+2_test.nc')
in_data[2] = o.to_array()

#points onto glacier
g_ind = np.where(in_data[4].values.flatten() > 0) 

now=time.time()
warnings.filterwarnings("ignore")

# loop over each point over the glacier
yearly_mb = pd.DataFrame(index=np.linspace(2011.68,2019.68,9))

  
#for i in g_ind[0]:
#    yearly_mb[f'i{i}']=0

# loop over each point over glacier
def parall(i):
    yearly_mb[f'i{i}']=0
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
            in_data[1].values.flatten()[i] # geometric altitude in y_omega+2, on glacier sfc
    
    melt_f = in_data[2].values.flatten()[i]
    # spinup
    pd_bucket0 = spinup(spin_years, altitude, melt_f)
    
    #get climate 2011.68-2020.68
    cl = get_climate_cal(years, altitude)
    
    pd_bucket = pd_bucket0
    
    a = 0
    suma = 0
    index=np.linspace(2011.68,2019.68,9)
    
    for yi in cl.index:
        m_b_tot = monthly_mb_sd(cl.loc[yi], pd_bucket)
        pd_bucket = m_b_tot[1]
        mbi = m_b_tot[0] 
        
        suma = mbi+suma # in mmwe
        if any(np.array(index)==yi): #if october
            yearly_mb[f'i{i}'][yi] = suma # yearly_sum
            #yearly_mb.append(suma)
            suma=0
    return yearly_mb

pool = mp.Pool(ncores) #4 cores

now=time.time()
yearly_mb = pool.map(parall, [i for i in g_ind[0]])
yearly_mb = pd.concat(yearly_mb, axis=1)
now1=time.time()
print(f'Time minimizing is: {abs(now - now1)}')

pool.close()

#initialize pd
#for i in g_ind[0]:
#    yearly_mb[f'i{i}']=0
#
#for i in g_ind[0]:
#    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
#    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
#            in_data[1].values.flatten()[i] # geometric altitude in y_omega+2, on glacier sfc
#    
#    melt_f = in_data[2].values.flatten()[i]
#    # spinup
#    pd_bucket0 = spinup(spin_years, altitude, melt_f)
#    
#    #get climate 2011.68-2020.68
#    cl = get_climate_cal(years, altitude)
#    
#    pd_bucket = pd_bucket0
#    
#    a = 0
#    suma = 0
#    index=np.linspace(2011.68,2019.68,9)
#    
#    for yi in cl.index:
#        m_b_tot = monthly_mb_sd(cl.loc[yi], pd_bucket)
#        pd_bucket = m_b_tot[1]
#        mbi = m_b_tot[0] 
#        
#        suma = mbi+suma # in mmwe
#        if any(np.array(index)==yi): #if october
#            yearly_mb[f'i{i}'][yi] = suma # yearly_sum
#            #yearly_mb.append(suma)
#            suma=0
#        
#    #yearly_mb[f'{i}'] =    
#    #residual.append(in_data[1].values.flatten()[i]-suma/rho)
#    print(f'observed dh mm: {in_data[1].values.flatten()[i]}, computed dh mm: {suma/rho}')
    
# plot
mass_b = []
for j in np.arange(0, 9):
    mass_b.append(yearly_mb.iloc[j].mean())

# geodetic mb:
glacier_mean = in_data[1].values.flatten()[g_ind[0]].mean()/9

plt.plot(np.linspace(2011.68,2019.68,9),mass_b,'go--' , label='computed_yearly_mb')
plt.hlines(np.mean(mass_b), 2011, 2020, colors='blue',label='computed_yearly_mb_mean')
plt.hlines(glacier_mean*1000*rho, 2011, 2020, colors='red',label='geodetic_mb')
plt.ylabel('yearly mass balance (mmwe)')
plt.xlabel('year')
#plt.xticks(np.linspace(2011.68,2019.68,9))
plt.legend()
plt.show()

   # print(f'Point altitude is: {altitude} and its melt year is: {melt_y_i}')

now1=time.time()
print(f'Time running future evolution is: {abs(now - now1)}')


#############################3
# ssp='ssp126'
# ensamble_name =  'ukesm1-0-ll_r1i1p1f2'
# a1=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}.nc')
# ensamble_name = 'gfdl-esm4_r1i1p1f1'

# a2=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}.nc')
# ensamble_name =  'ipsl-cm6a-lr_r1i1p1f1'

# a3=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}.nc')
# ensamble_name = 'mpi-esm1-2-hr_r1i1p1f1'

# a4=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}.nc')
# ensamble_name =  'mri-esm2-0_r1i1p1f1'
# a5=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_f_0_{y_alfa}-{y_omega}.nc')
# a_mean=(a1+a2+a3+a4+a5)/5

# a_mean=a_mean.to_array()

# # plot
# # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
# fig, ax1 = plt.subplots(figsize=(13, 3), ncols=1)
# # plot just the positive data and save the
# # color "mappable" object returned by ax1.imshow
# pos = ax1.imshow(a_mean.values[0], vmin=80, vmax=135)

# ax1.set_title(f'{cal}_{ssp}_mean_melt_f')
# fig.colorbar(pos, ax=ax1)
# fig.show()  
# ax1.set_ylabel('lat (pix)')
# ax1.set_xlabel('lon (pix)')
# plt.savefig(f'{out_path}/{n}x{n}/{cal}_{ssp}/mean_melt_factor_0_{y_alfa}-{y_omega}')
# plt.show()

# aaa = a_mean

# aaa.to_netcdf(f'{out_path}/{n}x{n}/{cal}_{ssp}/mean_melt_factor_0_{y_alfa}-{y_omega}.nc')

# ###########################
# #############################3
# ssp='ssp126'
# ensamble_name =  'ukesm1-0-ll_r1i1p1f2'
# a1=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_year_0_{y_alfa}-{y_omega}.nc')
# ensamble_name = 'gfdl-esm4_r1i1p1f1'

# a2=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_year_0_{y_alfa}-{y_omega}.nc')
# ensamble_name =  'ipsl-cm6a-lr_r1i1p1f1'

# a3=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_year_0_{y_alfa}-{y_omega}.nc')
# ensamble_name = 'mpi-esm1-2-hr_r1i1p1f1'

# a4=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_year_0_{y_alfa}-{y_omega}.nc')
# ensamble_name =  'mri-esm2-0_r1i1p1f1'
# a5=xr.open_dataset(f'{out_path}/{n}x{n}/{cal}_{ssp}/{ensamble_name}_melt_year_0_{y_alfa}-{y_omega}.nc')
# a_mean=(a1+a2+a3+a4+a5)/5

# a_mean=a_mean.to_array()

# # plot
# # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
# fig, ax1 = plt.subplots(figsize=(13, 3), ncols=1)
# # plot just the positive data and save the
# # color "mappable" object returned by ax1.imshow
# pos = ax1.imshow(a_mean.values[0][0], vmin=2020-1, vmax=2056)

# ax1.set_title(f'{cal}_{ssp}_melt_year_mean')
# fig.colorbar(pos, ax=ax1)
# ax1.set_ylabel('lat (pix)')
# ax1.set_xlabel('lon (pix)')
# fig.show()  
# plt.savefig(f'{out_path}/{n}x{n}/{cal}_{ssp}/mean_melt_year_0_{y_alfa}-{y_omega}')
# plt.show()

# aaa = a_mean

# aaa.to_netcdf(f'{out_path}/{n}x{n}/{cal}_{ssp}/mean_melt_year_0_{y_alfa}-{y_omega}.nc')