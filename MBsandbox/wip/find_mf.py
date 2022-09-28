#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calling functions from totalMB. To get a minimized MF
"""
from totalMB import (omnibus_minimize_mf, melt_f_update1)
from scipy.optimize import minimize_scalar
import numpy as np
from get_data_over_glacier import get_raster_data
import matplotlib.pyplot as plt

altitude = 3200  

res = minimize_scalar(omnibus_minimize_mf)
print(f'{int(res.x)} is the melt factor that minimizes abs(observed_mb-modelled_mb)')

result = []
for melt_f in np.arange(10,150,20):
    result.append(omnibus_minimize_mf(melt_f))

# TODO: do a minimiation keeping altitude cnstant and only haing melt_f as a variable
#omnibus_minimize_mf(melt_f, altitude)

#scipy.optimize.minimize(omnibus_minimize_mf(..., altitude=1000))

in_data = get_raster_data()
g_ind = np.where(in_data[4].values.flatten()>0) #points onto glacier

#altitude2011 = bed topo + thichness 2020 + thick diff


# loop over each point
m_f = []
for i in g_ind[0]:
    
    altitide = in_data[4].values.flatten()[i] + in_data[3].values.flatten()[i] + \
            in_data[1].values.flatten()[i]
    obs_mb = in_data[1].values.flatten()[i]
    
    obs_mb = obs_mb * 1000 #(in mm)

    res = minimize_scalar(omnibus_minimize_mf, args=(altitude, obs_mb)) # check mm and mmwe in side omnibus function
    m_f.append(res.x)

dum = np.zeros(100)
dum[g_ind] = m_f
dum_resh = np.reshape(dum, [10,10])

# fill melt_f in in_data list
in_data[2].values=dum_resh

# plot
plt.imshow(in_data[2])
plt.legend()
plt.show()

# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(20, 3), ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
pos = ax1.imshow(in_data[2], vmin=min(in_data[2].values.flatten()[in_data[2].values.flatten()>0]-10), vmax=max(in_data[2].values.flatten()))

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near
ax1.set_title('melt_factor, 2011-2020')
fig.colorbar(pos, ax=ax1)
