#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 11:08:42 2022

@author: francesc
"""

import xarray as xr
import params as params
from params import * 
#from projection_functions import omnibus_future
import pandas as pd
import params as params
import numpy as np
import matplotlib.pyplot as plt


#function:
def get_stats_proj(filepath):
    """
    Parameters
    ----------
    file : pd dataframe path, in form of pickle
        file from prediction. (.pkl)

    Returns
    -------
    rel_volume, rel_area, melt_year

    """
    obj = pd.read_pickle(filepath)
    area = [] 
    volume = []
    for i in obj.index: 
        area.append(np.float(sum(obj.loc[i]!=0)))
        volume.append(np.float(sum(obj.loc[i])))

    melt_y = [] 
    for j in obj.columns:
        try: 
            melt_y.append(2020.68+min(obj.index[obj[j]==0]))
        except:    
            print('some error message')
    
    
    return volume, area, melt_y


proj_rhos=[]
proj_enss=[]
proj_ssps=[]

for rho in [790, 850, 910]:
    for ensamble_name in ensamble_names[:-1]:
        for ssp in ssps:
            filepath=f'{out_path}/{n}x{n}/projection/df/{ensamble_name}_{ssp}_{int(rho)}_df.pkl'
            results=get_stats_proj(filepath)
            proj_ssps.append(results)
        proj_enss.append(proj_ssps)
        proj_ssps=[]
    proj_rhos.append(proj_enss)
    proj_enss=[]
proj_pd = pd.DataFrame(proj_rhos)       

obj=pd.read_pickle(f'{out_path}/{n}x{n}/projection/df/{ensamble_name}_{ssp}_{int(rho)}_df.pkl')

# plot

for i_ens in [0,1,2,3,4]:
    for i_ssp in [0,1,2]:
        if i_ssp==0 :
            ltype='b-'
        elif i_ssp==1 :
            ltype='g-'
        elif i_ssp==2 :
            ltype='r-'
        label={i_ens}
        plt.plot(np.array(obj.index[0:500]+2020.68),np.array(proj_pd.iloc[0][i_ens][i_ssp][1][0:500])/proj_pd.iloc[0][i_ens][i_ssp][1][0], ltype)
        plt.hlines(0.1, xmin=2020.68,xmax=2060)
        plt.title('relative area change 790')


for i_ens in [0,1,2,3,4]:
    for i_ssp in [0,1,2]:
        if i_ssp==0 :
            ltype='b-'
        elif i_ssp==1 :
            ltype='g-'
        elif i_ssp==2 :
            ltype='r-'
        label={i_ens}
        plt.plot(np.array(obj.index[0:500]+2020.68),np.array(proj_pd.iloc[0][i_ens][i_ssp][0][0:500])/proj_pd.iloc[0][i_ens][i_ssp][0][0], ltype)
        plt.hlines(0.1, xmin=2020.68,xmax=2060)
        plt.title('relative volume change 790')
  
#melt_y0=melt_y

# same for no sd results:
proj_rhos_nosd=[]
proj_enss_nosd=[]
proj_ssps_nosd=[]

for rho in [850]:
    for ensamble_name in ensamble_names[:-1]:
        for ssp in ssps:
            filepath=f'{out_path}/{n}x{n}/projection/df/{ensamble_name}_{ssp}_no_sd_df.pkl'
            results=get_stats_proj(filepath)
            proj_ssps_nosd.append(results)
        proj_enss_nosd.append(proj_ssps_nosd)
        proj_ssps_nosd=[]
    proj_rhos_nosd.append(proj_enss_nosd)
    proj_enss_nosd=[]
proj_pd_nosd = pd.DataFrame(proj_rhos_nosd)  


############################################
#ensemble_mean #rho=850
ensemble_sum=ensemble_sum850=ensemble_sum790=ensemble_sum910=ensemble_sum_nosd=0
for i_ens in [0,1,2,3,4]:
    for i_ssp in [0,1,2]:
        ensemble_sum790 = ensemble_sum790 + np.array(proj_pd.iloc[0][i_ens][i_ssp][2])
        ensemble_sum850 = ensemble_sum850 + np.array(proj_pd.iloc[1][i_ens][i_ssp][2])
        ensemble_sum910 = ensemble_sum910 + np.array(proj_pd.iloc[2][i_ens][i_ssp][2])
        ensemble_sum_nosd = ensemble_sum_nosd + np.array(proj_pd_nosd.iloc[0][i_ens][i_ssp][2])

ensemble_mean790=ensemble_sum790/15
ensemble_mean850=ensemble_sum850/15
ensemble_mean910=ensemble_sum910/15
ensemble_mean_nosd=ensemble_sum_nosd/15

ensemble_mean = ensemble_mean790-ensemble_mean850

#plot
g_ind = []
for s in obj.columns: g_ind.append(s[1:])

#plot 2D 
#proj_pd[0][0][0][2]

g_ind=np.array(g_ind, dtype=int)

dum = np.zeros(n * n)
dum[g_ind] = ensemble_mean_nosd-ensemble_mean850
dum_resh = np.reshape(dum, [n, n])

melt_y=dum_resh
melt_y[melt_y==0]=-1000

# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(13, 4), ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
pos = ax1.imshow(melt_y, vmin=-2, vmax=+20)

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

ax1.set_title('melt_year(rho=850kgm-3)-melt_year_nosd(rho=850kgm-3)')
plt.ylabel('latitude (m)')
plt.xlabel('longitude (m)')
fig.colorbar(pos, ax=ax1)
fig.show()


#n = params.n

#in_data = get_raster_data(n)

# open calibration m_f raster:
#o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_{int(rho)}.nc')
#in_data[2] = o.to_array()

#points onto glacier
#g_ind = np.where(in_data[4].values.flatten() > 0) 
    
#points onto glacier
g_ind = []
for s in obj.columns: g_ind.append(s[1:])

#plot 2D 
#proj_pd[0][0][0][2]

g_ind=np.array(g_ind, dtype=int)

for i_ens in [0,1,2,3,4]:
    for i_ssp in [0,1,2]:
        dum = np.zeros(n * n)
        dum[g_ind] = proj_pd.iloc[0][i_ens][i_ssp][2]
        dum_resh = np.reshape(dum, [n, n])

        melt_y=dum_resh
        
        # plot
        # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
        fig, ax1 = plt.subplots(figsize=(13, 3), ncols=1)
        # plot just the positive data and save the
        # color "mappable" object returned by ax1.imshow
        pos = ax1.imshow(melt_y, vmin=min(melt_y.flatten()[melt_y.flatten()>0]-5), vmax=max(melt_y.flatten()))

        # add the colorbar using the figure's method,
        # telling which mappable we're talking about and
        # which axes object it should be near

        ax1.set_title('melt_year')
        plt.ylabel('latitude (pix)')
        plt.xlabel('longitude (pix)')
        fig.colorbar(pos, ax=ax1)

#f'{out_path}/{n}x{n}/projection/{ensamble_name}_{ssp}_{int(rho*1000)}.pkl