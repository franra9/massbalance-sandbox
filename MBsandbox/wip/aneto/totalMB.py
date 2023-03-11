#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses the climate from 'climate.py' to compute the accumulated specific mass balance.
"""

import numpy as np
import pandas as pd
from params import rho, spin_years, wspinup

######### Compute melt factor using dummy data and surface distinction,
# without spinup period ########

# uptate melt factor according to surface type.
def melt_f_update1(time, melt_f_snow):
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
    accum - melt (mb, in mmwe, pd_bucket
    
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

    # age surface one month: # Warning: this ageing works only because the 71th position is also nan for our period and our glacier.
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
                
                melt_f = pd_non0buck['melt_factor'].iloc[0]
                dummelt = melt_f * temp4melt_m / fact

                # pd_bucket['mmwe'][ind] = - mbdiff = - prcpsol + meltf * t4m_lim

                temp4melt_limit = pd_non0buck['mmwe'].iloc[0] * fact / pd_non0buck['melt_factor'].iloc[0]
                
                if temp4melt_limit > temp4melt_m: # bucket at ind=ind not melted totally

                    pd_non0buck['mmwe'].iloc[0] = pd_non0buck['mmwe'].iloc[0] - dummelt # update bucket
                    
                    # Update melt
                    melt = melt + dummelt
                    
                    temp4melt_m = 0 # exit loop
                
                elif pd_non0buck['mmwe'].iloc[0] == 1.0: # melting oldest surface
  
                    melt = melt + dummelt
                    temp4melt_m = 0 # exit loop


                
                else:
                    # newest bucket melted
                    pd_non0buck.iloc[0] = np.nan
                    
                    # update bucket                    
                    pd_non0buck = pd_non0buck[np.isnan(pd_non0buck['mmwe']) == False]  
                    
                    # how much is melted in the youngest (now disappeared) bucket?
                    meltbuck = melt_f * temp4melt_limit / fact

                    # Update melt
                    melt = melt + meltbuck
                    
                    # Residual temp4melt
                    temp4melt_m = temp4melt_m - temp4melt_limit
           
    
    pd_bucket['mmwe'] = pd_non0buck['mmwe']
    
    #print([float(accum - melt), pd_bucket])
    #print(f'accum: {accum}')
    #print(f'melt: {melt}')
    return [float(accum - melt), pd_bucket] #in m/s # TODO: exit an updated pd_bucket

##############
# Create a function to put the melt_f as only parameter in order to 
# do the minimization afterwards
###############

from climate import get_climate, get_climate_cal

#def omnibus_minimize_mf(melt_f, cl=[0], altitude=0, obs_mb=0, years=np.array([2012, 2012.09])): #inicialize with some dummy values
#    """
#    This function is used to run amd minimize the melt factor to calibrate the sfc_mb model
#    
#    Parameters
#    ----------
#    melt_f: float 
#        melt factor from snow. Melt factor ice is then computed using the 
#        melt_f_update1 function.
#    
#    altitude: float 
#        altitude point
#    
#    obs_mb: float 
#        commited altitude change. Geometrical change (NOT swe)
#     
#    years : np array
#        float, natural float years of the period I want the climate.
#        
#    Returns
#    -------
#    Melt factor (float) that makes real ablation match observed ablation.
#
#    """
#    # initialize bucket:
#    pd_bucket = pd.DataFrame(index=np.linspace(0, 72, 73), columns = ['mmwe', 'melt_factor'])
#    pd_bucket['mmwe'] = np.nan
#    pd_bucket['mmwe'][-1] = 1
#    dum = []
#    
#    altitude=altitude
#    for i in pd_bucket.index:
#        dum.append(melt_f_update1(i, melt_f))  
#    # fill pd_bucket with melt factors
#    pd_bucket['melt_factor'] = dum
#
#    #initial values for climate # hardcoded for testing reasons
#    #years = np.round(np.linspace(2011, 2019, (2019-2011) * 12 + 1), 2) + 0.01 #add 0.01 to make sure there are not rounding errors
#
#    # apply monthly_mb_sd #
#    # get climate
#    #cl = get_climate(years, altitude)
#    cl = cl
#    
#    # monthly mb
#    total_m = []
#    summ_mb = 0
#    for iyr in np.arange(0, len(years)):
#        ## altitude change
#        # if summ_mb : altitude = altitude - mbbbbb (...)
#        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
#        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0]
#        total_m.append(mb)
#        summ_mb = summ_mb + mb
##    print(summ_mb - obs_mb)
##    print(melt_f)
##    return abs(summ_mb - obs_mb)
##    return summ_mb

def spinup(spin_years, altitude, melt_f):
    
    pd_bucket = pd.DataFrame(index=np.linspace(0, 72, 73), columns = ['mmwe', 'melt_factor'])
    pd_bucket['mmwe'] = np.nan
    pd_bucket['mmwe'][-1] = 1
    dum= []

    altitude=altitude
    for i in pd_bucket.index:
        dum.append(melt_f_update1(i, melt_f))  
    # fill pd_bucket with melt factors
    pd_bucket['melt_factor'] = dum
    
    cl = get_climate(spin_years, altitude)
    
    # monthly mb
    total_m = []
    summ_mb = 0
    summ_mb_i = 0
    
    for iyr in np.arange(0, len(spin_years)):
        ## altitude change
        #if abs(summ_mb_i) > 1000 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
        #    altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
        #    cl = get_climate(spin_years, altitude) # update climate for new altitude
        #    summ_mb_i = 0
        # if summ_mb : altitude = altitude - mbbbbb (...)
        
        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0]
        total_m.append(mb)
        summ_mb_i = summ_mb_i + mb
        summ_mb = summ_mb + mb

    return pd_bucket

def omnibus_minimize_mf(melt_f, altitude=0, obs_mb=0, years=np.array([2012, 2012.09])): #inicialize with some dummy values
    """
    Same as omnibus_minimize_mf but with altitude-climate change
    
    Parameters
    ----------
    melt_f: float 
        melt factor from snow. Melt factor ice is then computed using the 
        melt_f_update1 function.
       
    altitude: float 
        altitude point
    
    obs_mb: float 
        commited altitude change in mmwe.
     
    years : np array
        float, natural float years of the period I want the climate.
        
    Returns
    -------
    difference observation-computed mb
    """
    altitude0=altitude
    
    # spinup
    spin_years = np.round(np.linspace(2011-6, 2011, (2011 - (2011-6)) * 12 + 1), 2) + 0.68 #0.67 stands for october
    pd_bucket = spinup(spin_years, altitude, melt_f)
    
#    if wspinup == '':
#        pd_bucket = pd.DataFrame(index=np.linspace(0, 72, 73), columns = ['mmwe', 'melt_factor'])
#        pd_bucket['mmwe'] = np.nan
#        pd_bucket['mmwe'][-1] = 1
#        dum= []
#
#        altitude=altitude
#        for i in pd_bucket.index:
#            dum.append(melt_f_update1(i, melt_f))  
#        # fill pd_bucket with melt factors
#        pd_bucket['melt_factor'] = dum
        
    # start calibration, spinup done
    cl = get_climate_cal(years, altitude)

    # monthly mb
    total_m = []
    summ_mb = 0
    summ_mb_i = 0
    
    #reinicialize altitude to do calibration.
    altitude = altitude0
    
    for iyr in np.arange(0, len(years) + 12): # add 24 months to get to 2020.68
        ## altitude change
        #if abs(summ_mb_i/(1000 * rho)) > 1 : #5m w.e. # effect max ~4mm/month: (0.0298*4.5*30)
        #    altitude = altitude + (summ_mb_i / (1000 * rho))  #geometric m
        #    cl = get_climate_cal(years, altitude) # update climate for new altitude
        #    summ_mb_i = 0
            
        #update bucket
        pd_bucket = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[1]
        mb = monthly_mb_sd(cl.iloc[iyr], pd_bucket)[0]
        total_m.append(mb)
        summ_mb_i = summ_mb_i + mb
        summ_mb = summ_mb + mb
        
    if abs(summ_mb - obs_mb) < 100:
        print('-----------ini-------')
        print(f'g_ind: unknonw')
        print(f'residual:{summ_mb - obs_mb}')
        print(f'obs_mb:{obs_mb}')
        print(f'summ_mb:{summ_mb}')
        print(f'melt_factor:{melt_f}')
        print('----------end-------')
    return abs(summ_mb - obs_mb)

