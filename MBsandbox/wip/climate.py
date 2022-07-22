#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:45:17 2022
This script loads some climate that I can use for my purposes.
@author: francesc
"""
import numpy as np
import pandas as pd
import pointMB

years = np.round(np.linspace(2010, 2015, (2015-2010) * 12 + 1), 2)

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

# Test model for the priod 2010-2015  # years are in normal years, not hydro.
    for yr in years:
        clim.temp[yr], clim.temp4melt[yr], clim.prcp[yr], clim.prcpsol[yr] = \
            pointMB._get_climate(altitude, climate_type = 'monthly', year=yr)
    return clim

# clim is the climate in the 5 year test period
