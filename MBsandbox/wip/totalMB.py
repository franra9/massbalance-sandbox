#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses the climate from 'climate.py' to compute the accumulated specific mass balance.
"""
import numpy as np

def monthly_mb(climate, melt_f_bucket=None, melt_f=120):
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

    if melt_f_bucket is None: # temporary, first test.
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
    else:
        # bucket melt_f stuff
        print('TBA')
    return mb_month 

years = np.round(np.linspace(2010, 2015, (2015-2010) * 12 + 1), 2)
altitude = 3200

from climate import get_climate

cl = get_climate(years, altitude)

cl.loc[[2014]]

melt_f = 120

mb = 0

for yr in years:
    mb = mb + monthly_mb(cl.loc[yr])

from scipy.optimize import minimize_scalar
obs_mb = 5000 #mmweq
def mb_compare(melt_f):
    """
    compare obs and computed mb

    Parameters
    ----------
    melt_f : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

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




