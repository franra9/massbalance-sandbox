#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses the climate from 'climate.py' to compute the accumulated specific mass balance.
"""

import numpy as np
import pandas as pd

def monthly_mb(climate, melt_f):
    """
    Compute specific mass balance for a given climate.

    Parameters
    ----------
    climate : 
        t, temp2dformelt, prcp, prcpsol
    bucket : pandas dataframe
        accounts for the surface type
    
    Returns
    -------
    monthly mass balance

    """  
    # get 2D values, dependencies on height and time (days)
    t, temp4melt, prcp, prcpsol = climate.values
    # (days per month)
    # dom = 365.25/12  # len(prcpsol.T)
    fact = 12/365.25
    # attention, I should not use the days of years as the melt_f is
    # per month ~mean days of that year 12/daysofyear
    # to have the same unit of melt_f, which is
    # the monthly temperature sensitivity (kg /m² /mth /K),
    mb_daily = prcpsol - melt_f * temp4melt * fact

    mb_month = np.sum(mb_daily, axis=1)
    # more correct than using a mean value for days in a month
    #warnings.warn('there might be a problem with SEC_IN_MONTH'
    #              'as February changes amount of days inbetween the years'
    #              ' see test_monthly_glacier_massbalance()')
    return mb_month 

######### Compute melt factor using dummy data and no surface distinction ########
years = np.round(np.linspace(2010, 2015, (2015-2010) * 12 + 1), 2)
altitude = 3200

from climate import get_climate

cl = get_climate(years, altitude)
melt_f = 120
mb = 0

for yr in years:
    mb = mb + monthly_mb(cl.loc[yr], melt_f=melt_f)

from scipy.optimize import minimize_scalar
obs_mb = 7000 #mmweq


def mb_compare(melt_f):
    """
    compare obs and computed mb

    Parameters
    ----------
    melt_f : int
        melt factor.

    Returns
    -------
    float, difference between the computed and the observed mass balance

    """
    mb = 0
    years = np.round(np.linspace(2010, 2015, (2015-2010) * 12 + 1), 2)
    altitude = 3200
    cl = get_climate(years, altitude)

    for yr in years:    
        mb = mb + monthly_mb(cl.loc[yr], melt_f=melt_f)
    return abs(abs(obs_mb) + mb)

res = minimize_scalar(mb_compare)
print(f'{int(res.x)} is the melt factor that minimizes abs(observed_mb-modelled_mb)')

######### Compute melt factor using dummy data and surface distinction,
# without spinup period ########

# uptate melt factor according to surface type.
def melt_f_update(time):
    """
    Explonential decay of melt factor. 

    Parameters
    ----------
    time : int
        Age of the last nonzero bucket, in months

    Returns
    -------
    melt factor in kg /m2 /K /month

    """
    melt_f_snow = 100
    melt_f_ice = 2 * melt_f_snow # this is to have only one independent variable insteda of 2
    tau = 12 #tau_e_fold_yr=1 year as in Lilis code
    # add mmwe
    
    if time < 72:
        melt_f = melt_f_ice + (melt_f_snow - melt_f_ice) * np.exp(-time/tau) 
    else:
        melt_f = melt_f_ice
    
    return melt_f

def monthly_mb_sd(climate, pd_bucket):
    """
    Compute specific mass balance for a given climate, with surface distinction.

    Parameters
    ----------
    climate : 
        t, temp2dformelt, prcp, prcpsol
    bucket : pandas dataframe
        accounts for the surface type
    
    Returns
    -------
    monthly mass balance

    """  
    # get 2D values, dependencies on height and time (days)
    t, temp4melt, prcp, prcpsol = climate.values
    # (days per month)
    # dom = 365.25/12  # len(prcpsol.T)
    #fact = 12/365.25
    
    # check first non-null mmwe in pd_bucket
    temp4melt_m = sum(temp4melt[0])
    prcpsol_m = sum(prcpsol[0])
    
    # initialize melt and accum
    melt = 0
    accum = 0

    # age surface one month:
    aged_mmwe = np.roll(pd_bucket['mmwe'][0:71].values,1)
    pd_bucket['mmwe'][0:71] = aged_mmwe 
    
    # add solid ppt
    if prcpsol_m !=0:
        accum = prcpsol_m
        pd_bucket['mmwe'][:0] = prcpsol_m
            
 #  pd_non0buck = pd_bucket[pd_bucket['mmwe']!=0] or pd_bucket[np.isnan(pd_bucket['mmwe']) == False]
    pd_non0buck = pd_bucket[np.isnan(pd_bucket['mmwe']) == False]
        
    # now we want to find kow much do we melt.
        
    # melt term
    while temp4melt_m > 0: # melt bucket
            #if pd_non0buck[:0].index == 72: #ice bucket, the last one, no solid ppt
            # factor to corect month/day melt 
            fact = len(t[0])
                
            if len(pd_non0buck) == 1:    
                melt_f = pd_non0buck['melt_factor']
                melt = melt + melt_f * temp4melt_m / fact
                              
                temp4melt_m = 0
                
            else: # go to the youngest bucket
                
                melt_f = pd_non0buck['melt_factor'][0]
                dummelt = melt_f * temp4melt_m / fact

                # pd_bucket['mmwe'][ind] = - mbdiff = - prcpsol + meltf * t4m_lim

                temp4melt_limit = pd_bucket['mmwe'][0] * fact / pd_bucket['melt_factor'][0]
                
                if temp4melt_limit > temp4melt_m: # bucket at ind=ind not melted totally
                    melt = melt + dummelt
                    pd_non0buck['mmwe'][0] = pd_non0buck['mmwe'][0] - melt # update bucket
                    
                    temp4melt_m = 0 # exit loop
                
                else:
                    # newest bucket melted
                    pd_non0buck.iloc[0] = np.nan
                    
                    # update bucket                    
                    pd_non0buck = pd_non0buck[np.isnan(pd_non0buck['mmwe']) == False]                    
                    meltbuck = melt_f * temp4melt_limit / fact

                    #pd_non0buck = pd_non0buck[1:] # the equivalent thing is done three lines above
                    melt = melt + meltbuck
                    
                    temp4melt_m = temp4melt_m - temp4melt_limit
           
    
    pd_bucket['mmwe'] = pd_non0buck['mmwe']
    
    print([float(accum - melt), pd_bucket])
    print(f'accum: {accum}')
    print(f'melt: {melt}')
    return [float(accum - melt), pd_bucket] #in m/s # TODO: exit an updated pd_bucket

    # attention, I should not use the days of years as the melt_f is
    # per month ~mean days of that year 12/daysofyear
    # to have the same unit of melt_f, which is
    # the monthly temperature sensitivity (kg /m² /mth /K),
    #mb_daily = prcpsol - melt_f * temp4melt * fact

    #mb_month = np.sum(mb_daily, axis=1)
    # more correct than using a mean value for days in a month
    #warnings.warn('there might be a problem with SEC_IN_MONTH'
    #              'as February changes amount of days inbetween the years'
    #              ' see test_monthly_glacier_massbalance()')
    #return mb_month 

# initialize bucket:
pd_bucket = pd.DataFrame(index=np.linspace(0, 72, 73), columns = ['mmwe', 'melt_factor'])
pd_bucket['mmwe'] = np.nan
pd_bucket['mmwe'][-1] = 1
dum= []

for i in pd_bucket.index:
    dum.append(melt_f_update(i))  
# fill pd_bucket with melt factors
pd_bucket['melt_factor'] = dum

#initial values for climate
years = np.round(np.linspace(2010, 2015, (2015-2010) * 12 + 1), 2) + 0.01 #add 0.01 to make sure there are not rounding errors
altitude = 3200

from climate import get_climate

cl = get_climate(years, altitude)
melt_f = 120
mb = 0

for yr in years:
    mb = mb + monthly_mb(cl.loc[yr], melt_f=melt_f)

from scipy.optimize import minimize_scalar
obs_mb = 7000 #mmweq

# apply monthly_mb_sd #
# get climate
cl = get_climate(years, altitude)

#pd_bucket = monthly_mb_sd(cl.iloc[0], pd_bucket)[1]

######## test initialization for 1 year:
pd_bucket = monthly_mb_sd(cl.iloc[0], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[1], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[2], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[3], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[4], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[5], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[6], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[7], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[8], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[9], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[10], pd_bucket)[1]
pd_bucket = monthly_mb_sd(cl.iloc[11], pd_bucket)[1]




