#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 2.

Perform calibration checks using w5e5 and ensemble data.
calibration period results for all ssp scenarios

@author: francesc
"""
import xarray as xr
import params as params
from get_data_over_glacier import get_raster_data
from projection_functions import omnibus_future
from params import *
import pandas as pd
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
o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_{int(rho*1000)}.nc')
in_data[2] = o.to_array()

#points onto glacier
g_ind = np.where(in_data[4].values.flatten() > 0) 

now=time.time()
warnings.filterwarnings("ignore")

yearly_mb = pd.DataFrame(index=np.linspace(0,8,9))

# loop over each point over the glacier
for i in g_ind[0]:
    yearly_mb[f'i{i}']=0

# loop over each point over glacier
########################
def parall(i):
#for i in g_ind[0][:]:
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')

    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
             in_data[1].values.flatten()[i] # geometric altitude in y_omega+2, on glacier sfc
    obs_mb = in_data[1].values.flatten()[i]
    melt_f = in_data[2].values.flatten()[i]
    
    altitude0=altitude

    # spinup
    spin_years = np.round(np.linspace(2011-6, 2011, (2011 - (2011-6)) * 12 + 1), 2) + 0.68 #0.67 stands for october
    pd_bucket = spinup(spin_years, altitude, melt_f)
    
        
    # start calibration, spinup done
    cl = get_climate_cal(years, altitude)

    # monthly mb
    total_m = []
    summ_mb = 0
    summ_mb_i = 0
    summ_mb_ii = 0
    #reinicialize altitude to do calibration.
    altitude = altitude0
    
    for iyr in np.arange(0, len(years) + 12): # add 24 months to get to 2020.68
        ## altitude change
        if abs(summ_mb_i) > 5000 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
            altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
            cl = get_climate_cal(years, altitude) # update climate for new altitude
            summ_mb_i = 0
        
        #update bucket
        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0]
        total_m.append(mb)
        summ_mb_ii = summ_mb_ii + mb
        summ_mb_i = summ_mb_i + mb
        summ_mb = summ_mb + mb
#        if (iyr+1) % 12 == 0:
#            print(summ_mb_ii)
#            print(iyr)
#            yearly_mb[f'i{i}'][-1+(iyr+1)/12] = summ_mb_ii # yearly_sum
#            #yearly_mb.append(summ_mb)
#            summ_mb_ii = 0
    #yearly_mb[f"i{i}"] = np.array([total_m]).reshape(9,12).sum(axis=1)
    
    print('-----------ini-------')
    print(f'g_ind: {i}')
    print(f'residual:{summ_mb - obs_mb*1000*0.85}')
    print(f'obs_mb:{obs_mb*1000*0.85}')
    print(f'summ_mb:{summ_mb}')
    print(f'melt_factor:{melt_f}')
    print(f'summ years:{yearly_mb[f"i{i}"].sum()}')
    print('----------end-------')
    
    return total_m

pool = mp.Pool(ncores) #4 cores

now=time.time()
out = pool.map(parall, [i for i in g_ind[0][:]])
##yearly_mb = pd.concat(yearly_mb, axis=1)

now1=time.time()
print(f'Time minimizing is: {abs(now - now1)}')

pool.close()

for ii,i in enumerate(g_ind[0]):
    yearly_mb[f'i{i}']=np.array([out[ii]]).reshape(9,12).sum(axis=1)

# plot
mass_b = []
for j in np.arange(0, 9):
    mass_b.append(yearly_mb.iloc[j].mean())

# geodetic mb:
glacier_mean = in_data[1].values.flatten()[g_ind[0]].mean()/9

#plt.plot(mass_b,'go--' , label='computed_yearly_mb')
plt.plot(np.linspace(2012,2020,9),mass_b,'go--' , label='computed_yearly_mb')
plt.hlines(np.mean(mass_b), 2012, 2020, colors='blue',label='computed_yearly_mb_mean')
plt.hlines(glacier_mean*1000*rho, 2012, 2020, colors='red',label='geodetic_mb')
plt.ylabel('yearly mass balance (mmwe)')
plt.xlabel('ydro year')
#plt.xticks(np.linspace(2011.68,2019.68,9))
plt.legend()
plt.savefig(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb')
#plt.show()

yearly_mb.to_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb.pkl')
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

