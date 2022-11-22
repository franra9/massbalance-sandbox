#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 11:58:52 2022

@author: francesc

Script to plot different ssp scenarios and w5e5 data
"""
import pandas as pd
import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt
data_path = '/home/francesc/data/aneto_glacier/climate/OGGM/OGGM-sfc-type/per_glacier/RGI60-11/RGI60-11.03/RGI60-11.03208/'

s00=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_ukesm1-0-ll_r1i1p1f2_ssp126_no_correction.nc')).to_pandas()
s01=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_ukesm1-0-ll_r1i1p1f2_ssp370_no_correction.nc')).to_pandas()
s02=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_ukesm1-0-ll_r1i1p1f2_ssp585_no_correction.nc')).to_pandas()
s10=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mri-esm2-0_r1i1p1f1_ssp126_no_correction.nc')).to_pandas()
s11=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mri-esm2-0_r1i1p1f1_ssp370_no_correction.nc')).to_pandas()
s12=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mri-esm2-0_r1i1p1f1_ssp585_no_correction.nc')).to_pandas()
s20=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_gfdl-esm4_r1i1p1f1_ssp126_no_correction.nc')).to_pandas()
s21=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_gfdl-esm4_r1i1p1f1_ssp370_no_correction.nc')).to_pandas()
s22=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_gfdl-esm4_r1i1p1f1_ssp585_no_correction.nc')).to_pandas()
s30=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_ipsl-cm6a-lr_r1i1p1f1_ssp126_no_correction.nc')).to_pandas()
s31=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_ipsl-cm6a-lr_r1i1p1f1_ssp370_no_correction.nc')).to_pandas()
s32=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_ipsl-cm6a-lr_r1i1p1f1_ssp585_no_correction.nc')).to_pandas()
s40=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mpi-esm1-2-hr_r1i1p1f1_ssp126_no_correction.nc')).to_pandas()
s41=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mpi-esm1-2-hr_r1i1p1f1_ssp370_no_correction.nc')).to_pandas()
s42=xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mpi-esm1-2-hr_r1i1p1f1_ssp585_no_correction.nc')).to_pandas()
w5e5=xr.open_dataset(os.path.join(data_path, 'climate_historical_daily_W5E5.nc')).to_pandas()

ind=(s00.index<'2020-10-01') * (s00.index>'2011-09-30')

s00=s00[ind]
s01=s01[ind]
s02=s02[ind]
s10=s10[ind]
s11=s11[ind]
s12=s12[ind]
s20=s20[ind]
s21=s21[ind]
s22=s22[ind]
s30=s30[ind]
s31=s31[ind]
s32=s32[ind]
s40=s40[ind]
s41=s41[ind]
s42=s42[ind]
w5e5=w5e5[(w5e5.index<'2020-10-01') * (w5e5.index>'2011-09-30')]

s00=s00.resample('y').mean().temp
s01=s01.resample('y').mean().temp
s02=s02.resample('y').mean().temp
s10=s10.resample('y').mean().temp
s11=s11.resample('y').mean().temp
s12=s12.resample('y').mean().temp
s20=s20.resample('y').mean().temp
s21=s21.resample('y').mean().temp
s22=s22.resample('y').mean().temp
s30=s30.resample('y').mean().temp
s31=s31.resample('y').mean().temp
s32=s32.resample('y').mean().temp
s40=s40.resample('y').mean().temp
s41=s41.resample('y').mean().temp
s42=s42.resample('y').mean().temp
w5e5=w5e5.resample('y').mean().temp

xr.open_dataset(os.path.join(data_path, 'gcm_data_daily_ISIMIP3b_mri-esm2-0_r1i1p1f1_ssp126_no_correction.nc'))

plt.plot(w5e5,'go--' , label=f'calib')
plt.plot(s00,'r-' , label=f'126')
plt.plot(s01,'g-' , label=f'370')
plt.plot(s02,'b-' , label=f'585')
plt.plot(s10,'r-' )
plt.plot(s11,'g-' )
plt.plot(s12,'b-')
plt.plot(s20,'r-' )
plt.plot(s21,'g-')
plt.plot(s22,'b-' )
plt.plot(s30,'r-')
plt.plot(s31,'g-' )
plt.plot(s32,'b-')
plt.plot(s40,'r-')
plt.plot(s41,'g-' )
plt.plot(s42,'b-')

#plt.hlines(calib_mean, 2011, 2020, color='black',label='cal_mean')
#plt.hlines(glacier_mean*1000*rho, 2011, 2020, colors='brown',label='geodetic_mb')
plt.ylabel('mean yearly temp')
plt.xlabel('year')
#plt.xticks(np.linspace(2011.68,2019.68,9))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.savefig(f'{out_path}/{n}x{n}/calibration/cliamte_comp_test')
plt.show()
