#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:46:30 2022
prediction of total melt. Exits a netcdf and png with the melt year for each point.
@author: francesc
"""
import xarray as xr
import params as params
from params import * 
from get_data_over_glacier import get_raster_data
#from projection_functions import omnibus_future
import pandas as pd
import params as params
from totalMB import (omnibus_minimize_mf, spinup,monthly_mb_sd)
from scipy.optimize import minimize_scalar
import numpy as np
from get_data_over_glacier import get_raster_data
import matplotlib.pyplot as plt
import time
import warnings
import multiprocessing as mp
from climate import get_climate

# get data from rasters
n = params.n

in_data = get_raster_data(n)

# open calibration m_f raster:
#o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_{int(rho*1000)}.nc')
o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_no_sd.nc')
#w5e5_melt_f_0_2011-2019+1_no_sd.nc
in_data[2] = o.to_array()

#points onto glacier
g_ind = np.where(in_data[4].values.flatten() > 0) 

now=time.time()
warnings.filterwarnings("ignore")

yearly_mb = pd.DataFrame(index=np.linspace(0,79,80))

# loop over each point over the glacier
for i in g_ind[0]:
    yearly_mb[f'i{i}']=0

# loop over each point over glacier
########################
def parallll(i):
#for i in g_ind[0][:]:
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')

    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
             in_data[1].values.flatten()[i] # geometric altitude in y_omega, on glacier sfc
    obs_mb = in_data[1].values.flatten()[i]
    melt_f = in_data[2].values.flatten()[i]
    thick20 = in_data[3].values.flatten()[i]
    
    altitude0=altitude

    # spinup
    spin_years = np.round(np.linspace(2020-6, 2020, (2020 - (2020-6)) * 12 + 1), 2) + 0.68 #0.67 stands for october
    pd_bucket = spinup(spin_years, altitude, melt_f)
    
        
    # start calibration, spinup done
    cl = get_climate(years_fcst, altitude)

    # monthly mb
    total_ice = [thick20]
    summ_mb = 0
    summ_mb_i = 0
    #reinicialize altitude to do calibration.
    altitude = altitude0
    
    for iyr in np.arange(0, len(years_fcst)): # add 24 months to get to 2020.68
        ## altitude change
        if abs(summ_mb_i/(1000*rho)) > 5 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
            altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
            cl = get_climate(years_fcst, altitude) # update climate for new altitude
            summ_mb_i = 0
        
        #update bucket
        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0]
        summ_mb_i = summ_mb_i + mb
        summ_mb = summ_mb + mb

        total_ice.append(thick20 + (summ_mb / (1000 * rho)))
        #print(total_ice[-1])        
        #print(thick20 + total_mm / (1000 * rho) )
        #print(thick20 + total_mm / (1000 * rho) < 0)
        if thick20 + summ_mb / (1000 * rho) < 0: #mm + mmwe/(1000rho)=mm+mm
            dum=np.zeros(80*12)
            dum[:len(total_ice)] = total_ice
            break
    
    print('-----------ini-------')
    print(f'g_ind: {i}')
    #print(f'residual:{summ_mb - obs_mb*1000*0.85}')
    #print(f'obs_mb:{obs_mb*1000*0.85}')
    #print(f'summ_mb:{summ_mb}')
    #print(f'melt_factor:{melt_f}')
    #print(f'summ years:{yearly_mb[f"i{i}"].sum()}')
    print(f'melt year: {2020+iyr/12}')
    print('----------end-------')
    
    return dum

pool = mp.Pool(ncores) #4 cores

now=time.time()
out = pool.map(parallll, [i for i in g_ind[0][:]])
##yearly_mb = pd.concat(yearly_mb, axis=1)

now1=time.time()
print(f'Time minimizing is: {abs(now - now1)}')

pool.close()

monthly_ev = pd.DataFrame(index=np.linspace(0,79,80*12))

for ii,i in enumerate(g_ind[0][:]):
    monthly_ev[f'i{i}']=np.array([out[ii]])[0]#.reshape(80,12)#.sum(axis=1)

#sanity_check:
monthly_ev=np.clip(monthly_ev,0,50)

remaining_ice = monthly_ev.sum(axis=1)

#sanity_check:
remaining_ice=np.clip(remaining_ice, 0, 1999999)

thick20_0 = in_data[3].values.flatten()[g_ind[0]]
thick20 = np.clip(thick20_0, 0, 199999)

#2020 is thichness 2020
#remaining_ice = np.roll(remaining_ice,1)
remaining_ice[0] = thick20.sum()


plt.plot(remaining_ice[:]/thick20.sum(), label=f'{ensamble_name}_{ssp}')
plt.xlabel('years from 2020')
plt.ylabel('Normalized volume [-]')
plt.legend()
#plt.savefig(f'{out_path}/{n}x{n}/projection/{ensamble_name}_{ssp}_{int(rho*1000)}')
plt.savefig(f'{out_path}/{n}x{n}/projection/{ensamble_name}_{ssp}_no_sd')

evol=remaining_ice[:]/thick20.sum()
#evol.to_pickle(f'{out_path}/{n}x{n}/projection/{ensamble_name}_{ssp}_{int(rho*1000)}.pkl')
evol.to_pickle(f'{out_path}/{n}x{n}/projection/{ensamble_name}_{ssp}_no_sd.pkl')

now1=time.time()
print(f'Time running future evolution is: {abs(now - now1)}')

#save data frame
#monthly_ev.to_pickle(f'{out_path}/{n}x{n}/projection/df/{ensamble_name}_{ssp}_{int(rho*1000)}_df.pkl')
monthly_ev.to_pickle(f'{out_path}/{n}x{n}/projection/df/{ensamble_name}_{ssp}_no_sd_df.pkl')


#d=xr.DataArray(out)
#d.to_pickle(f'{out_path}/{n}x{n}/projection/{ensamble_name}_{ssp}_df.pkl')

# # get data from rasters
# n = params.n

# in_data = get_raster_data(n)

# # open calibration m_f raster:
# o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1.nc')
# in_data[2] = o.to_array()

# #points onto glacier
# g_ind = np.where(in_data[4].values.flatten() > 0) 

# now=time.time()
# warnings.filterwarnings("ignore")

# # loop over each point over the glacier
# melt_y = []
# for i in g_ind[0]:
    
#     print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
#     altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
#             in_data[1].values.flatten()[i] # geometric altitude in y_omega, on glacier sfc
    
#     melt_y_i = omnibus_future(melt_f = in_data[2].values.flatten()[i], altitude = altitude, \
#                    thick20 = in_data[3].values.flatten()[i], years_fcst = years_fcst)
#     melt_y.append(melt_y_i)
#     print(f'Point altitude is: {altitude} and its melt year is: {melt_y_i}')

#now1=time.time()
#print(f'Time running future evolution is: {abs(now - now1)}')
#
## grid size
plt_map=False
if plt_map:
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
    pos = ax1.imshow(in_data[2].values[0], vmin=min(in_data[2].values[0].flatten()[in_data[2].values[0].flatten()>0]-1), vmax=max(in_data[2].values[0].flatten()))
    
    # add the colorbar using the figure's method,
    # telling which mappable we're talking about and
    # which axes object it should be near
    ax1.set_title(f'{y_alfa}-{y_omega}_{n}x{n}_{cal}_{ssp}_{ensamble_name}')
    fig.colorbar(pos, ax=ax1)
    fig.show()  
    plt.savefig(f'{out_path}/{n}x{n}/w5e5_/{ensamble_name}_{ssp}_melt_year_0_{y_alfa}-{y_omega}')
    plt.show()
    
    aaa = in_data[2]
    
    aaa.to_netcdf(f'{out_path}/{n}x{n}/w5e5_/{ensamble_name}_{ssp}_melt_year_0_{y_alfa}-{y_omega}.nc')

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