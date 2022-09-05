#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calling functions from totalMB. To get a minimized MF
"""
from totalMB import (omnibus_minimize_mf, melt_f_update1)
from scipy.optimize import minimize_scalar
import numpy as np

altitude = 3200  

res = minimize_scalar(omnibus_minimize_mf)
print(f'{int(res.x)} is the melt factor that minimizes abs(observed_mb-modelled_mb)')


altitude = 3400  

res = minimize_scalar(omnibus_minimize_mf)
print(f'{int(res.x)} is the melt factor that minimizes abs(observed_mb-modelled_mb)')
    
result = []
for melt_f in np.arange(10,150,20):
    result.append(omnibus_minimize_mf(melt_f))

# TODO: do a minimiation keeping altitude cnstant and only haing melt_f as a variable
#omnibus_minimize_mf(melt_f, altitude)

#scipy.optimize.minimize(omnibus_minimize_mf(..., altitude=1000))