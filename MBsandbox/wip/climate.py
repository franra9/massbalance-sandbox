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

years = np.linspace(2010, 2015, (2015-2010) * 12 + 1) 

column_names = ["temp", "temp4melt", "prcp", "prcpsol"]

df = pd.DataFrame(columns = column_names, index=years)

# Test model for the priod 2010-2015  # years are in normal years, not hydro.
year=np.array([2010.]) 
temp0, tempformelt,prcp,prcsol = pointMB._get_climate(3000, climate_type = 'monthly', year=year)

year=np.array([2011.0])
temp1, tempformelt,prcp,prcsol = pointMB._get_climate(3000, climate_type = 'monthly', year=year)

year=np.array([2011.09])
temp2, tempformelt,prcp,prcsol = pointMB._get_climate(3000, climate_type = 'monthly', year=year)

year=np.array([2011.18])
temp3, tempformelt,prcp,prcsol = pointMB._get_climate(3000, climate_type = 'monthly', year=year)

df.temp[2010] = temp0