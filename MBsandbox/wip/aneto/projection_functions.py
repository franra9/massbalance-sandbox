#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create functions 

@author: franra9 11102022
"""
import xarray as xr
import params as params
from params import * 
from get_data_over_glacier import get_raster_data
import pandas as pd
from totalMB import *


#in_data = get_raster_data(n)
# open calibration m_f raster:
#in_data[2] = xr.open_dataset(f'{out_path}/{y_alfa}-{y_omega}_{n}x{n}_{cal}_{ssp}_melt_f_0.nc')
spin_years = params.spin_years

#def spinup(spin_years, altitude):
#    
#    pd_bucket = pd.DataFrame(index=np.linspace(0, 72, 73), columns = ['mmwe', 'melt_factor'])
#    pd_bucket['mmwe'] = np.nan
#    pd_bucket['mmwe'][-1] = 1
#    dum= []
#
#    altitude=altitude
#    for i in pd_bucket.index:
#        dum.append(melt_f_update1(i, melt_f))  
#    # fill pd_bucket with melt factors
#    pd_bucket['melt_factor'] = dum
#    
#    cl = get_climate(spin_years, altitude)
#    
#    # monthly mb
#    total_m = []
#    summ_mb = 0
#    summ_mb_i = 0
#    
#    for iyr in np.arange(0, len(spin_years)):
#        ## altitude change
#        if abs(summ_mb_i) > 5000 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
#            altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
#            cl = get_climate(years, altitude) # update climate for new altitude
#            summ_mb_i = 0
#        # if summ_mb : altitude = altitude - mbbbbb (...)
#        
#        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
#        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0]
#        total_m.append(mb)
#        summ_mb_i = summ_mb_i + mb
#        summ_mb = summ_mb + mb
#
#    return pd_bucket


def omnibus_future(melt_f, altitude=0, thick20=1, years_fcst=np.array([2012, 2012.09])): #inicialize with some dummy values
    """
    This function is used to forecast point glacier melt
    
    Parameters
    ----------
    melt_f : float 
        melt factor from snow, coming from calibration 092011-092020. 
    
    pd_bucket0: pd dataframe
        initial conditions of the ice bucket.
    
    altitude : float 
        altitude point, in m
    
    thick20 : float
        Ice thickness 2020 from GPR`, in m, (geometric)
    
    years_fcst : np array
        float, natural float years of the period I want to forecast the glacier evolution.
        
    Returns
    -------
    Float year of single point melt.

    """
    years = years_fcst
    altitude0=altitude
    if wspinup == 'wspinup':
        spin_years = np.linspace(years[0]-6, years[0], 72)
        pd_bucket = spinup(spin_years, altitude, melt_f)
    
    pd_bucket = pd.DataFrame(index=np.linspace(0, 72, 73), columns = ['mmwe', 'melt_factor'])
    pd_bucket['mmwe'] = np.nan
    pd_bucket['mmwe'][-1] = 1
    dum= []

    altitude=altitude
    for i in pd_bucket.index:
        dum.append(melt_f_update1(i, melt_f))  
        # fill pd_bucket with melt factors
        pd_bucket['melt_factor'] = dum
        
    # apply monthly_mb_sd #
    # get climate
    cl = get_climate(years, altitude)
    
    # monthly mb
    total_m = []
    total_mm = 0
    summ_mb = 0 # counter for each 5mwe melt
    mb = 0
    
    for iyr in np.arange(0, len(years)):
        ## altitude change
        if abs(summ_mb/1000*rho) > 5 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
            altitude = altitude + (summ_mb / (1000 * rho))  #geometric m
            cl = get_climate(years, altitude) # update climate for new altitude
            summ_mb = 0
            #print('climate updated')
            
        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0] # in mmwe
        #print(mb)
        #total_m.append(mb)
        total_mm = total_mm + mb
        summ_mb = summ_mb + mb
        
        #print(thick20 + total_mm / (1000 * rho) )
        #print(thick20 + total_mm / (1000 * rho) < 0)
        if thick20 + total_mm / (1000 * rho) < 0: #mm + mmwe/(1000rho)=mm+mm
            break


    #print(altitude)
    return years[iyr]
