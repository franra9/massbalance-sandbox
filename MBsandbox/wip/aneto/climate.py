#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:45:17 2022
This script loads some climate that I can use for my purposes.
@author: francesc
"""
import pandas as pd
import pointMB
import numpy as np
#import aneto_projections_bayescalibration

# calibration
#years = np.round(np.linspace(2011, 2019, (2019-2011) * 12 + 1), 2)

# prediction
#years = np.round(np.linspace(2011, 2013, (2013-2011) * 12 + 1), 2)


######## get climate ########
def get_climate(years, altitude):
    """
    gets daily data from the float years in 'year', organized monthly.

    Parameters
    ----------
    year : np array
        float, natural years of the period I want the climate
    
    altitude : int
        altitude m.a.s.l. where to get the climate
        

    Returns
    -------
    pd dataframe: 2m? temperature at my altitude, temperature4melt, t
    otal? precipitation, solid precipitation per day, grouped by month.

    """
    
    column_names = ["temp", "temp4melt", "prcp", "prcpsol"]

    clim = pd.DataFrame(columns = column_names, index=years)

# Test model for the priod 2011-2019  # years are in normal years, not hydro.
    for yr in years:
        clim.temp[yr], clim.temp4melt[yr], clim.prcp[yr], clim.prcpsol[yr] = \
            pointMB._get_climate(altitude, climate_type = 'monthly', year=yr)
    return clim

######## get climate ########
def get_climate_cal(years, altitude):
    """
    gets daily data from the float years in 'year', organized monthly.
    extra year is added to meet glaciological calibration period (2011.68-2020.68)

    Parameters
    ----------
    year : np array
        float, natural years of the period I want the climate (2011.68-2018.68)
    
    altitude : int
        altitude m.a.s.l. where to get the climate
        

    Returns
    -------
    pd dataframe: 2m? temperature at my altitude, temperature4melt, t
    otal? precipitation, solid precipitation per day, grouped by month.

    """
    
    column_names = ["temp", "temp4melt", "prcp", "prcpsol"]

    clim = pd.DataFrame(columns = column_names, index=years)

# Test model for the priod 2011-2019  # years are in normal years, not hydro.
    for yr in years:
        clim.temp[yr], clim.temp4melt[yr], clim.prcp[yr], clim.prcpsol[yr] = \
            pointMB._get_climate(altitude, climate_type = 'monthly', year=yr)
    
    new_years = np.append(years,np.round(np.linspace(2019.76, 2020.68,12),2))
    new_clim = pd.DataFrame(columns = column_names, index=new_years)
    new_clim.iloc[-12:] = clim.iloc[-12:]
    new_clim.iloc[:96] = clim
    return new_clim
