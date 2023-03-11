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
o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1.nc')#{int(rho*1000)}.nc')
in_data[2] = o.to_array()

#points onto glacier
g_ind = np.where(in_data[4].values.flatten() > 0) 

now=time.time()
warnings.filterwarnings("ignore")

# loop over each point over the glacier
yearly_mb = pd.DataFrame(index=np.linspace(0,8,9))
  

for i in g_ind[0]:
    yearly_mb[f'i{i}']=0

#for i in g_ind[0][:]:
def paralll(i):
    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
             in_data[1].values.flatten()[i] # geometric altitude in y_omega+2, on glacier sfc
    obs_mb = in_data[1].values.flatten()[i]
    melt_f = in_data[2].values.flatten()[i]
    # spinup
    spin_years = np.round(np.linspace(2011-6, 2011, (2011 - (2011-6)) * 12 + 1), 2) + 0.68 #0.67 stands for october
    pd_bucket = spinup(spin_years, altitude, melt_f)
    
    altitude0=altitude

       
    # start calibration, spinup done
    cl = get_climate(years, altitude)
    
    # monthly mb
    total_m = []
    summ_mb = 0
    summ_mb_i = 0
    summ_mb_ii = 0
    #reinicialize altitude to do calibration.
    altitude = altitude0
    
    for iyr in np.arange(0, len(years)): # add 24 months to get to 2020.68
        ## altitude change
        if abs(summ_mb_i/(1000*rho)) > 5 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
            altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
            cl = get_climate(years, altitude) # update climate for new altitude
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
# loop over each point over glacier
#def parall(i):
#    yearly_mb[f'i{i}']=0
#    print(f'{np.where(g_ind[0]==i)} out of {len(g_ind[0])} points')
#    altitude = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
#            in_data[1].values.flatten()[i] # geometric altitude in y_omega+2, on glacier sfc
#    
#    melt_f = in_data[2].values.flatten()[i]
#    # spinup
#    pd_bucket0 = spinup(spin_years, altitude, melt_f)
#    
#    #get climate 2011.68-2020.68
#    cl = get_climate(years, altitude)
#    
#    pd_bucket = pd_bucket0
#    
#    tot = 0
#    suma = 0
#    index=np.linspace(2011.68,2019.68,9)
#    summ_mb_i = 0
#    for yi in cl.index:
#        ## altitude change
#        if abs(summ_mb_i) > 5000 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
#            altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
#            cl = get_climate(years, altitude) # update climate for new altitude
#            summ_mb_i = 0
#        
#        mb = monthly_mb_sd(cl.loc[yi], pd_bucket)[0]
#        pd_bucket = monthly_mb_sd(cl.loc[yi], pd_bucket)[1]
#        
#        summ_mb_i = summ_mb_i + mb
#        # add surface altitude change!
#        
#        suma = mb+suma # in mmwe
#        if any(np.array(index)==yi): #if october
#            yearly_mb[f'i{i}'][yi] = suma # yearly_sum
#            #yearly_mb.append(suma)
#            suma=0
#    return yearly_mb

pool = mp.Pool(ncores) #4 cores

now=time.time()
out = pool.map(paralll, [i for i in g_ind[0][:]])
now1=time.time()
print(f'Time minimizing is: {abs(now - now1)}')

pool.close()

for ii,i in enumerate(g_ind[0]):
    yearly_mb[f'i{i}']=np.array([out[ii]]).reshape(9,12).sum(axis=1)

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

plt.plot(np.linspace(2011.68,2019.68,9),mass_b,'go--' , label=f'{ensamble_name}_{ssp}_yr_mb')
#plt.hlines(np.mean(mass_b), 2011, 2020, colors='blue',label='comp_yr_mb_mean')
plt.hlines(glacier_mean*1000*rho, 2011, 2020, colors='red',label='geodetic_mb')
plt.ylabel('yearly mass balance (mmwe)')
plt.xlabel('year')
#plt.xticks(np.linspace(2011.68,2019.68,9))
plt.legend()
plt.savefig(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}')
#plt.show()

yearly_mb.to_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
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
#######################

ensamble_name=ensamble_names[0]

ssp=ssps[0]

r00=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r00=r00.mean(axis=1)

ensamble_name=ensamble_names[0]

ssp=ssps[1]

r01=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r01=r01.mean(axis=1)
ensamble_name=ensamble_names[0]

ssp=ssps[2]

r02=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r02=r02.mean(axis=1)
ensamble_name=ensamble_names[1]

ssp=ssps[0]

r10=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r10=r10.mean(axis=1)
ensamble_name=ensamble_names[1]

ssp=ssps[1]

r11=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r11=r11.mean(axis=1)
ensamble_name=ensamble_names[1]

ssp=ssps[2]

r12=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r12=r12.mean(axis=1)
ensamble_name=ensamble_names[2]

ssp=ssps[0]

r20=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r20=r20.mean(axis=1)
ensamble_name=ensamble_names[2]

ssp=ssps[1]

r21=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r21=r21.mean(axis=1)
ensamble_name=ensamble_names[2]

ssp=ssps[2]

r22=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r22=r22.mean(axis=1)
ensamble_name=ensamble_names[3]

ssp=ssps[0]

r30=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r30=r30.mean(axis=1)
ensamble_name=ensamble_names[3]

ssp=ssps[1]

r31=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r31=r31.mean(axis=1)
ensamble_name=ensamble_names[3]

ssp=ssps[2]

r32=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r32=r32.mean(axis=1)
ensamble_name=ensamble_names[4]

ssp=ssps[0]

r40=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r40=r40.mean(axis=1)
ensamble_name=ensamble_names[4]

ssp=ssps[1]

r41=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r41=r41.mean(axis=1)
ensamble_name=ensamble_names[4]

ssp=ssps[2]

r42=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r42=r42.mean(axis=1)

calib=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/w5e5__yearly_mb.pkl')
calib=calib.mean(axis=1)
calib_mean=calib.mean()

#plot:
    
plt.plot(np.linspace(2012.68,2020.68,9),calib,'go--' , label=f'calib')
plt.plot(np.linspace(2012.68,2020.68,9),r00,'r-' , label=f'126')
plt.plot(np.linspace(2012.68,2020.68,9),r01,'g-' , label=f'370')
plt.plot(np.linspace(2012.68,2020.68,9),r02,'b-' , label=f'585')
plt.plot(np.linspace(2012.68,2020.68,9),r10,'r-' )
plt.plot(np.linspace(2012.68,2020.68,9),r11,'g-' )
plt.plot(np.linspace(2012.68,2020.68,9),r12,'b-')
plt.plot(np.linspace(2012.68,2020.68,9),r20,'r-' )
plt.plot(np.linspace(2012.68,2020.68,9),r21,'g-')
plt.plot(np.linspace(2012.68,2020.68,9),r22,'b-' )
plt.plot(np.linspace(2012.68,2020.68,9),r30,'r-')
plt.plot(np.linspace(2012.68,2020.68,9),r31,'g-' )
plt.plot(np.linspace(2012.68,2020.68,9),r32,'b-')
plt.plot(np.linspace(2012.68,2020.68,9),r40,'r-')
plt.plot(np.linspace(2012.68,2020.68,9),r41,'g-' )
plt.plot(np.linspace(2012.68,2020.68,9),r42,'b-')

plt.hlines(calib_mean, 2012.68, 2020.68, color='black',label='geodeticMB')#label='cal_mean')
#plt.hlines(glacier_mean*1000*rho, 2011.68, 2020.68, colors='brown',label='geodetic_mb')
plt.ylabel('yearly mass balance (mmwe)')
plt.xlabel('year')
#plt.xticks(np.linspace(2011.68,2019.68,9))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f'{out_path}/{n}x{n}/calibration/all_yearly_mb_{int(rho*1000)}')
#plt.show()

#yearly_mb.to_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb.pkl')
    # print(f'Point altitude is: {altitude} and its melt year is: {melt_y_i}')

now1=time.time()
print(f'Time running future evolution is: {abs(now - now1)}')