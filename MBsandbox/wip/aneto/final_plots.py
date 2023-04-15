#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for the final version:
@froura 03012022
"""

import xarray as xr
import params as params
from params import * 
#from projection_functions import omnibus_future
import pandas as pd
import params as params
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

###############################################
##########  plots calibration  ################
###############################################

#cmap #https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
#cmap = plt.cm.jet  # define the colormap
# extract all colors from the .jet map
#cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#cmaplist[0] = (.5, .5, .5, 1.0)

# create the new map
#cmap = mpl.colors.LinearSegmentedColormap.from_list(
#    'Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
#bounds = np.linspace(0, 20, 21)
#norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

####
# previous stuff:
from get_data_over_glacier import get_raster_data
n = params.n

in_data = get_raster_data(n)
g_ind = np.where(in_data[4].values.flatten() > 0) #points onto glacier

# read files
rmf_790_5 = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_790.nc')
rmf_850_5 = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_850.nc')
rmf_910_5 = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_910.nc')

#rmf_850_5_t = xr.open_dataset(f'{out_path}/{n}x{n}/calibration_/w5e5_melt_f_0_2011-2019+1.nc')#H=None


rmf_850_0 = xr.open_dataset(f'{out_path}/{n}x{n}/calibration_/w5e5_melt_f_0_2011-2019+1.nc')
rmf_850_1 = xr.open_dataset(f'{out_path}/{n}x{n}/calibration_1m/w5e5_melt_f_0_2011-2019+1_1m.nc')

rmf_850_5_no_sd = xr.open_dataset(f'{out_path}/{n}x{n}/calibration_/w5e5_melt_f_0_2011-2019+1_no_sd.nc')

rmf_850_0_new = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_2011-2019+1_850_none_new.nc')
rmf_850_1_new = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_2011-2019+1_850_1m_new.nc')

#put data in an easier variable, onto the glacier already
dum = np.zeros(n * n)
dum[g_ind] = rmf_790_5['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_790_5 = dum_resh

# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_850_5['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_850_5 = dum_resh

# grid size
#dum = np.zeros(n * n)
#dum[g_ind] = rmf_850_5_t['__xarray_dataarray_variable__'].values.flatten()[g_ind]
#dum_resh = np.reshape(dum, [n, n])
#mf_850_5_t = dum_resh

# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_910_5['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_910_5 = dum_resh


# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_850_1['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_850_1 = dum_resh

# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_850_1_new['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_850_1_new = dum_resh


# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_850_0['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_850_0 = dum_resh

# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_850_0_new['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_850_0_new = dum_resh

# grid size
dum = np.zeros(n * n)
dum[g_ind] = rmf_850_5_no_sd['__xarray_dataarray_variable__'].values.flatten()[g_ind]
dum_resh = np.reshape(dum, [n, n])

mf_850_5_no_sd = dum_resh


#####################
# preliminary plots #
#####################
vmin=1.6
vmax=2
plt.imshow(mf_850_5, vmin=vmin,vmax=vmax); plt.title('H=5m'); plt.colorbar()
plt.show()
plt.imshow(mf_850_1, vmin=vmin,vmax=vmax); plt.title('H=1m'); plt.colorbar()
plt.show()
plt.imshow(mf_850_1_new, vmin=vmin,vmax=vmax); plt.title('H=1m,new'); plt.colorbar()
plt.show()
plt.imshow(mf_850_0, vmin=vmin,vmax=vmax); plt.title('H=None'); plt.colorbar()
plt.show()
plt.imshow(mf_850_0_new, vmin=vmin,vmax=vmax); plt.title('H=None,new'); plt.colorbar()
plt.show()
plt.imshow(mf_850_5_no_sd, vmin=vmin,vmax=vmax); plt.title('H=5m, nosurfdist'); plt.colorbar()
plt.show()


##################################
# preliminary plots - comparison #
##################################
vmax = 0.2
vmin = -0.6
plt.imshow((mf_850_1-mf_850_0), vmin=vmin,vmax=vmax); plt.title('H=1m-H=None'); plt.colorbar()
plt.show()
plt.imshow((mf_850_5-mf_850_0), vmin=vmin,vmax=vmax); plt.title('H=5m-H=None'); plt.colorbar()
plt.show()
plt.imshow((mf_850_5-mf_850_1), vmin=vmin,vmax=vmax); plt.title('H=5m-H=1m'); plt.colorbar()
plt.show()
plt.imshow((mf_850_0-mf_850_5), vmin=vmin,vmax=vmax); plt.title('H=None-H=5m'); plt.colorbar()
plt.show()
plt.imshow((mf_850_1-mf_850_5), vmin=vmin,vmax=vmax); plt.title('H=1m-H=5m'); plt.colorbar()
plt.show()

vmax = 0.8
vmin = -0.2
plt.imshow((mf_850_1-mf_850_0)/(mf_850_0), vmin=vmin,vmax=vmax); plt.title('H=1m-H=None %'); plt.colorbar()
plt.show()
plt.imshow((mf_850_0-mf_850_5)/(mf_850_5), vmin=vmin,vmax=vmax); plt.title('H=5m-H=None %'); plt.colorbar()
plt.show()
plt.imshow((mf_850_1-mf_850_5)/(mf_850_5), vmin=vmin,vmax=vmax); plt.title('H=1m-H=5m %'); plt.colorbar()
plt.show()

plt.imshow(mf_850_0/30.4); plt.title('H=None'); plt.colorbar()
plt.show()
plt.imshow(mf_850_5_no_sd/30.4); plt.title('H=5m, nosurfdist'); plt.colorbar()
plt.show()
plt.imshow(mf_850_5/30.4); plt.title('H=5m'); plt.colorbar()
plt.show()
plt.imshow(mf_850_1/30.4); plt.title('H=1m'); plt.colorbar()
plt.show()
plt.imshow(mf_850_0/30.4); plt.title('H=None'); plt.colorbar()
plt.show()

plt.imshow(mf_850_0_new/30.4); plt.title('H=None,new'); plt.colorbar()
plt.show()
plt.imshow(mf_850_5_no_sd/30.4); plt.title('H=5m, nosurfdist'); plt.colorbar()
plt.show()
plt.imshow(mf_850_5/30.4); plt.title('H=5m'); plt.colorbar()
plt.show()
plt.imshow(mf_850_1_new/30.4); plt.title('H=1m,new'); plt.colorbar()
plt.show()
plt.imshow(mf_850_0_new/30.4); plt.title('H=None,new'); plt.colorbar()
plt.show()


# get latlon for the plot axes
x0 = rmf_790_5.x.values
y0 = rmf_790_5.y.values

#get latlon with origin in my glacier
x = x0 - x0[0]
y = y0 - y0[-1]
default_x_ticks = range(len(x))
default_y_ticks = range(len(y))
x = x[np.round(np.linspace(0,39,5),0).astype(int)].astype(int)
x[-1] = x[-1]+0.5
y = y[np.round(np.linspace(0,39,5),0).astype(int)].astype(int)
y[0] = y[0]+0.5
default_x_ticks = np.array(default_x_ticks)
default_x_ticks = default_x_ticks[np.round(np.linspace(0,39,5),0).astype(int)]
default_y_ticks = np.array(default_y_ticks)
default_y_ticks = default_y_ticks[np.round(np.linspace(0,39,5),0).astype(int)]


###########
#Figure 3.1: Melt factor field over Aneto glacier, in mm(day)−1K−1 in a 32mx27m
#resolution using ρice = 850kgm−3 and H = 5m
toplot = mf_850_5/30.4
toplot[toplot==0] = np.nan
#remove 3 spots
toplot[np.isclose(toplot,0)]=np.nan

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define boundaries
vmax = 2
vmin = 1.3#min(toplot.flatten()[toplot.flatten()>0])

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,15)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor for ' + r'$\rho_{ice}$'+'=850kgm'+r'$^{-3}$'+', H=5m '+r'$[mmd^{-1}K^{-1}$]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor ($mmd^{-1}K^{-1}$)")

plt.savefig(f'{out_path}/{n}x{n}/plots/fig31.png')
plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########
#Figure 3.2: μ field change for H = 5m. μ(ρice = 790kgm−3) − μ(ρice = 850kgm−3) (a),
#and μ(ρice = 910kgm−3 − ρice = 850kgm−3) (b), in mm(d)−1K−1.   
toplot_a = mf_790_5 - mf_850_5
toplot_b = mf_910_5 - mf_850_5
toplot_a = toplot_a/30.4
toplot_b = toplot_b/30.4


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define boundaries
vmax = 0.05
vmin = -0.05 #min(toplot.flatten()[toplot.flatten()>0])

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.PiYG
bounds = np.linspace(vmin,vmax,21)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_a, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor change\n' + r'$\mu(\rho_{ice}=790kgm^{-3})-\mu(\rho_{ice}=850kgm^{-3})$' + r'$[mmd^{-1}K^{-1}$]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change ($mmd^{-1}K^{-1}$)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig32a.png')
plt.show()

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.PiYG
bounds = np.linspace(vmin,vmax,21)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_b, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor change\n' + r'$\mu(\rho_{ice}=910kgm^{-3})-\mu(\rho_{ice}=850kgm^{-3})$' + r'$[mmd^{-1}K^{-1}$]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change ($mmd^{-1}K^{-1}$)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig32b.png')
plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########
#Figure 3.3: Relative change in % of the temperature sensitivity field over Aneto glacier
#after calibrating with ρice = 790kgm3 instead of ρice = 850kgm3, at
#32mx27m resolution.
toplot = 100 * (mf_790_5 - mf_850_5)/mf_850_5 # ULL!!! problem at 3 points
points_outliers_i = toplot < -4
toplot[toplot < -4] = np.nan

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define boundaries
vmax = 1
vmin = -3#min(toplot.flatten()[toplot.flatten()>0])

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,17)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('relative melt factor change\n' + r'$\mu(\rho_{ice}=790kgm^{-3})-\mu(\rho_{ice}=850kgm^{-3})$' + ' [%]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change (%)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig33.png')

plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Figure 3.4: μ field absolute change μ(Hec = None) − μ(Hec = 15) and μ(Hec =
#1) − μ(Hec = 5m) (b), In mm(day)−1K−1 for ρice = 850 at 32mx27m resolution.
#toplot_a = mf_850_1 - mf_850_0
#toplot_b = mf_850_5 - mf_850_0
toplot_a = mf_850_0 - mf_850_5
toplot_b = mf_850_1 - mf_850_5

toplot_a[toplot_a == 0] = np.nan
toplot_b[toplot_b == 0] = np.nan
toplot_a=toplot_a/30.4
toplot_b=toplot_b/30.4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define boundaries
#vmax = -0.00
#vmin = -0.05 #min(toplot.flatten()[toplot.flatten()>0])
vmax = 0.06
vmin = 0.00 #min(toplot.flatten()[toplot.flatten()>0])

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,13)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_a, norm=norm, cmap=cmap)

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor change\n' + r'$\mu(H=None)-\mu(H=5m)$' + r' [$mmd^{-1}K^{-1}$]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)

plt.savefig(f'{out_path}/{n}x{n}/plots/fig34a.png')
plt.show()

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,13)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_b, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor change\n' + r'$\mu(H=1m)-\mu(H=5m)$' + r' [$mmd^{-1}K^{-1}$]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change ($mmd^{-1}K^{-1}$)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig34b.png')
plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Figure 3.5: μ field relative change μ(Hec = N one) − μ(Hec = 1m) (a) and μ(Hec =
#none) − μ(Hec = 5m) (b), In % for ρice = 850 at 32mx27m resolution.
toplot_a = 100 * (mf_850_0 - mf_850_5)/mf_850_5 #PROBLEMA AQUI
toplot_b = 100 * (mf_850_1 - mf_850_5)/mf_850_5 #PROBLEMA AQUI
toplot_b[points_outliers_i] = np.nan
toplot_a[points_outliers_i] = np.nan

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define boundaries
vmax = 3
vmin = 0 #min(toplot.flatten()[toplot.flatten()>0])

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,16)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_a, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor relative change\n' + r'$\mu(H=None)-\mu(H=5m)$' + r' [%]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change (%)")

plt.savefig(f'{out_path}/{n}x{n}/plots/fig35a.png')
plt.show()

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,16)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_b, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor relative change\n' + r'$\mu(H=1m)-\mu(H=5m)$' + r' [%]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change (%)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig35b.png')
plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Figure 3.6: μ field absolute change (a), in mmd−1K−1 and relative change (b) in %,
#omitting surface distinction and H = 5m.
toplot_a = (mf_850_5_no_sd - mf_850_5)/30.4
toplot_b = 100 * (mf_850_5_no_sd - mf_850_5)/mf_850_5

toplot_a[toplot_a==0] = np.nan
toplot_a[toplot_a>-0.20] = np.nan


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define boundaries
vmax = -0.26
vmin = -0.38 #min(toplot.flatten()[toplot.flatten()>0])
#vmax=-0.3
#vmin=-30

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,16)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_a, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor change\n' + r'$\mu(without)-\mu(with)$' + ' surface distinction' + r' [$mmd^{-1}K^{-1}$]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor ($mmd^{-1}K^{-1}$)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig36a.png')
plt.show()

# plot #
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
vmax = -15
vmin = -24 #min(toplot.flatten()[toplot.flatten()>0])


fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)#figsize=(7, 7), ncols=1)
# Use imshow to plot the data
#pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,21)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot_b, norm=norm, cmap=cmap)#cmap='PiYG')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title('melt factor change\n' + r'$\mu(without)-\mu(with)$' + ' surface distinction' + r' [%]')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')
#plt.colorbar(pos, ax=ax1)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1)#, label="melt factor change (%)")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig36b.png')
plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Figure 3.7: Correlation between altitude and temperature sensitivity (ρice = 850kgm3,
#H = 5m).
sf_altitude_ind = in_data[3].values.flatten()+in_data[4].values.flatten()
sf_altitude_flat = sf_altitude_ind[sf_altitude_ind>1000]

toplot = mf_850_5.flatten()[sf_altitude_ind>1000]/30.4

plt.scatter(sf_altitude_flat, toplot, marker='+')
plt.title('melt factor as a function of altitude')
plt.ylabel('melt factor ' + r'$[mmd^{-1}K^{-1}]$')
plt.xlabel('altitude [masl]')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(sf_altitude_flat, toplot, 1)
#add linear regression line to scatterplot 
plt.plot(sf_altitude_flat, m*sf_altitude_flat+b, c="r")
plt.savefig(f'{out_path}/{n}x{n}/plots/fig37.png')
plt.show()
#TODO: add regression line!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############3 ULLLLLLLLLLLLLLLLLL change n a 10 hardcoded ##############

n=10

#Figure 3.8: Yearly SMB corresponding to the calibration climate dataset (w5e5) and the
#15 climate forcing combinations coming from ISIMIP3b (5 GCMs, 3 ssp each) averaged
#over the whole glacier.

ensamble_name=ensamble_names[0]

ssp=ssps[0]

r00=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r00=r00.mean(axis=1)

ensamble_name=ensamble_names[0]

ssp=ssps[1]

r01=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r01=r01.mean(axis=1)
ensamble_name=ensamble_names[0]

ssp=ssps[2]

r02=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r02=r02.mean(axis=1)
ensamble_name=ensamble_names[1]

ssp=ssps[0]

r10=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r10=r10.mean(axis=1)
ensamble_name=ensamble_names[1]

ssp=ssps[1]

r11=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r11=r11.mean(axis=1)
ensamble_name=ensamble_names[1]

ssp=ssps[2]

r12=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r12=r12.mean(axis=1)
ensamble_name=ensamble_names[2]

ssp=ssps[0]

r20=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r20=r20.mean(axis=1)
ensamble_name=ensamble_names[2]

ssp=ssps[1]

r21=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r21=r21.mean(axis=1)
ensamble_name=ensamble_names[2]

ssp=ssps[2]

r22=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r22=r22.mean(axis=1)
ensamble_name=ensamble_names[3]

ssp=ssps[0]

r30=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r30=r30.mean(axis=1)
ensamble_name=ensamble_names[3]

ssp=ssps[1]

r31=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r31=r31.mean(axis=1)
ensamble_name=ensamble_names[3]

ssp=ssps[2]

r32=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r32=r32.mean(axis=1)
ensamble_name=ensamble_names[4]

ssp=ssps[0]

r40=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r40=r40.mean(axis=1)
ensamble_name=ensamble_names[4]

ssp=ssps[1]

r41=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r41=r41.mean(axis=1)
ensamble_name=ensamble_names[4]

ssp=ssps[2]

r42=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb_{int(rho*1000)}.pkl')
r42=r42.mean(axis=1)

calib=pd.read_pickle(f'{out_path}/{n}x{n}/calibration/w5e5__yearly_mb.pkl')
calib=calib.mean(axis=1)
calib_mean=calib.mean()

#plot:
    
plt.plot(np.linspace(2012.68,2020.68,9),r00,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r01,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r02,'0.8' , label=f'ISIMIP members')
plt.plot(np.linspace(2012.68,2020.68,9),r10,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r11,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r12,'0.8')
plt.plot(np.linspace(2012.68,2020.68,9),r20,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r21,'0.8')
plt.plot(np.linspace(2012.68,2020.68,9),r22,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r30,'0.8')
plt.plot(np.linspace(2012.68,2020.68,9),r31,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r32,'0.8')
plt.plot(np.linspace(2012.68,2020.68,9),r40,'0.8')
plt.plot(np.linspace(2012.68,2020.68,9),r41,'0.8' )
plt.plot(np.linspace(2012.68,2020.68,9),r42,'0.8')

plt.hlines(calib_mean, 2012.68, 2020.68, color='black',label='geodeticMB')#label='cal_mean')

isimip_mean=(r00.sum() + r01.sum()+r02.sum()+r10.sum()+r11.sum()+r12.sum()+r20.sum()+r21.sum()+
                    r22.sum()+r30.sum()+r31.sum()+r32.sum()+r40.sum()+r41.sum()+r42.sum())/(15*9)
plt.hlines(isimip_mean, 2012.68, 2020.68, color='red',label='ensemble mean')#label='cal_mean')
#plt.hlines(glacier_mean*1000*rho, 2011.68, 2020.68, colors='brown',label='geodetic_mb')
plt.ylabel('yearly mass balance [mmwe]')
plt.xlabel('calendar year')
#plt.xticks(np.linspace(2011.68,2019.68,9))
plt.plot(np.linspace(2012.68,2020.68,9),calib,'go--' , label=f'calibration')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f'{out_path}/{10}x{10}/calibration/all_yearly_mb_{int(rho*1000)}')
plt.savefig(f'{out_path}/40x40/plots/fig38.png')
plt.show()

#yearly_mb.to_pickle(f'{out_path}/{n}x{n}/calibration/{ensamble_name}_{ssp}_yearly_mb.pkl')
# print(f'Point altitude is: {altitude} and its melt year is: {melt_y_i}')

#now1=time.time()
#print(f'Time running future evolution is: {abs(now - now1)}')

######################################################
################ plots prediction ####################
######################################################
import xarray as xr
import time
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

#ensemble_mean = ensemble_mean790-ensemble_mean850

#plot
# g_ind = []
# for s in obj.columns: g_ind.append(s[1:])

# #plot 2D 
# #proj_pd[0][0][0][2]

# g_ind=np.array(g_ind, dtype=int)

# dum = np.zeros(n * n)
# dum[g_ind] = ensemble_mean_nosd-ensemble_mean850
# dum_resh = np.reshape(dum, [n, n])

# melt_y=dum_resh
# melt_y[melt_y==0]=-1000

# # plot
# # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
# fig, ax1 = plt.subplots(figsize=(13, 4), ncols=1)
# # plot just the positive data and save the
# # color "mappable" object returned by ax1.imshow
# pos = ax1.imshow(melt_y, vmin=-2, vmax=+20)

# # add the colorbar using the figure's method,
# # telling which mappable we're talking about and
# # which axes object it should be near

# ax1.set_title('melt_year(rho=850kgm-3)-melt_year_nosd(rho=850kgm-3)')
# ax1.set_aspect(1090.5/1273.5)
# plt.ylabel('latitude from origin latitude [m]')
# plt.xlabel('longitude from origin longitude [m]')

# fig.colorbar(pos, ax=ax1)
# fig.show()


# #n = params.n

# #in_data = get_raster_data(n)

# # open calibration m_f raster:
# #o = xr.open_dataset(f'{out_path}/{n}x{n}/w5e5_/w5e5_melt_f_0_{y_alfa}-2019+1_{int(rho)}.nc')
# #in_data[2] = o.to_array()

# #points onto glacier
# #g_ind = np.where(in_data[4].values.flatten() > 0) 
    
# #points onto glacier
# g_ind = []
# for s in obj.columns: g_ind.append(s[1:])

# #plot 2D 
# #proj_pd[0][0][0][2]

# g_ind=np.array(g_ind, dtype=int)

# for i_ens in [0,1,2,3,4]: #ensemble member
#     for i_ssp in [0,1,2]: #ssp
#         dum = np.zeros(n * n)
#         dum[g_ind] = proj_pd.iloc[0][i_ens][i_ssp][2]
#         dum_resh = np.reshape(dum, [n, n])

#         melt_y=dum_resh
        
#         # plot
#         # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
#         fig, ax1 = plt.subplots(figsize=(13, 3), ncols=1)
#         # plot just the positive data and save the
#         # color "mappable" object returned by ax1.imshow
#         pos = ax1.imshow(melt_y, vmin=min(melt_y.flatten()[melt_y.flatten()>0]-5), vmax=max(melt_y.flatten()))

#         # add the colorbar using the figure's method,
#         # telling which mappable we're talking about and
#         # which axes object it should be near

#         ax1.set_title('melt_year')
#         ax1.set_aspect(1090.5/1273.5)
#         plt.ylabel('latitude from origin latitude [m]')
#         plt.xlabel('longitude from origin longitude [m]')

#         fig.colorbar(pos, ax=ax1)


################################
######Template
################################
# #define boundaries
# vmax = 2
# vmin = 1.3#min(toplot.flatten()[toplot.flatten()>0])

# # plot #
# # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
# fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# # Use imshow to plot the data
# #pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
# cmap = mpl.cm.viridis
# bounds = np.linspace(vmin,vmax,15)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

# pos = ax1.imshow(toplot, norm=norm, cmap=cmap)#cmap='PiYG')

# plt.xticks(default_x_ticks, x, rotation=-45)
# plt.yticks(default_y_ticks, y)
# ax1.set_title('melt factor for ' + r'$\rho_{ice}$'+'=850kgm'+r'$^{-3}$'+', H=5m '+r'$[mmd^{-1}K^{-1}$]')
# ax1.set_aspect(1090.5/1273.5)
# plt.ylabel('latitude from origin latitude [m]')
# plt.xlabel('longitude from origin longitude [m]')
# #plt.colorbar(pos, ax=ax1)

# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              ax=ax1)#, label="melt factor ($mmd^{-1}K^{-1}$)")

# plt.savefig(f'{out_path}/{n}x{n}/plots/fig31.png')
# plt.show()

######################################
######################################
######################################

#Figure 3.9: Melt year for each point independently, with ρice = 850kgm3, H=5m. En-
#semble mean for all GCMs and ssp scenarios

# vmax = 2
# vmin = 1.3#min(toplot.flatten()[toplot.flatten()>0])

# # plot #
# # https://matplotlib.org/stable/gallery/color/colorbar_basics.html
# fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# # Use imshow to plot the data
# #pos = ax1.imshow(toplot, vmin=min(toplot.flatten()[toplot.flatten()>0]), vmax=2, cmap=cmap)#cmap='PiYG')
# cmap = mpl.cm.viridis
# bounds = np.linspace(vmin,vmax,15)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

# pos = ax1.imshow(toplot, norm=norm, cmap=cmap)#cmap='PiYG')

# plt.xticks(default_x_ticks, x, rotation=-45)
# plt.yticks(default_y_ticks, y)
# ax1.set_title('melt factor for ' + r'$\rho_{ice}$'+'=850kgm'+r'$^{-3}$'+', H=5m '+r'$[mmd^{-1}K^{-1}$]')
# ax1.set_aspect(1090.5/1273.5)
# plt.ylabel('latitude from origin latitude [m]')
# plt.xlabel('longitude from origin longitude [m]')
# #plt.colorbar(pos, ax=ax1)

# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              ax=ax1)#, label="melt factor ($mmd^{-1}K^{-1}$)")

# plt.savefig(f'{out_path}/{n}x{n}/plots/fig31.png')
# plt.show()


# Figure 3.9: Melt year for each point independently, with ρice = 850kgm3, H=5m. En-
# semble mean for all GCMs and ssp scenarios.

dum = np.zeros(n * n)
dum[g_ind] = ensemble_mean850 #[0], [ensemble], [ssp], [790, 850, 910]
toplot = np.reshape(dum, [n, n])

tmp_limits = toplot < 2020
toplot[toplot < 2020] = np.nan

toplot[points_outliers_i] = np.nan

melt_y=toplot 
        
# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
pos = ax1.imshow(toplot, vmin=2020, vmax=2055)





vmin=2020
vmax=2055
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,2056,10)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot, norm=norm, cmap=cmap)
# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

ax1.set_title(f'melt year for ' + r'$\rho=850kgm^{-3} [yr]$')
ax1.set_aspect(1090.5/1273.5)
plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)

#fig.colorbar(pos, ax=ax1)


plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1)#norm=norm, cmap=cmap),ax=ax1)
plt.savefig(f'{out_path}/40x40/plots/fig39.png')
plt.show()

#Figure 3.10: Melt year for each point independently, with ρice = 790kgm3 (a), ρice =
#910kgm3 (b), H=5m.

# TODO: canviar colors, com el que he fet servir abans per diferenciar <0 i >0

#a)

dum = np.zeros(n * n)
dum[g_ind] = ensemble_mean790 - ensemble_mean850# proj_pd.iloc[0][i_ens][i_ssp][2]
toplot = np.reshape(dum, [n, n])

toplot[tmp_limits] = np.nan

toplot[points_outliers_i] = np.nan

melt_y=toplot
        
# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
#pos = ax1.imshow(melt_y, vmin=-2, vmax=+2)
vmin=-2
vmax=+2
cmap = cmap = mpl.cm.PiYG
bounds = np.linspace(vmin,vmax,17)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot, norm=norm, cmap='PiYG')

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title(f'melt year difference [yr]\n' + r'$meltYear(\rho=790kgm^{-3})-meltYear(\rho=850kgm^{-3})$')
ax1.set_aspect(1090.5/1273.5)

plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')

#fig.colorbar(pos, ax=ax1)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1)#norm=norm, cmap=cmap),ax=ax1)
plt.savefig(f'{out_path}/40x40/plots/fig310a.png')
plt.show()


# b)
#dum = np.zeros(n * n)
#dum[g_ind] = ensemble_mean910 - ensemble_mean850
#toplot = np.reshape(dum, [n, n])
#
#toplot[tmp_limits] = np.nan
#
#toplot[points_outliers_i] = np.nan
#
#melt_y=toplot
#        
## plot
## https://matplotlib.org/stable/gallery/color/colorbar_basics.html
#fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
## plot just the positive data and save the
## color "mappable" object returned by ax1.imshow
#pos = ax1.imshow(melt_y, vmin=-2, vmax=+2)
#
# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near
#
#ax1.set_title(f'melt year difference\n' + r'$meltYear(\rho=910kgm^{-3})-meltYear(\rho=850kgm^{-3})$')
#ax1.set_aspect(1090.5/1273.5)
#plt.ylabel('latitude from origin latitude [m]')
#plt.xlabel('longitude from origin longitude [m]')
#
#plt.xticks(default_x_ticks, x, rotation=-45)
#plt.yticks(default_y_ticks, y)
#
#fig.colorbar(pos, ax=ax1)
#plt.show()

dum = np.zeros(n * n)
dum[g_ind] = ensemble_mean910 - ensemble_mean850
toplot = np.reshape(dum, [n, n])

toplot[tmp_limits] = np.nan

toplot[points_outliers_i] = np.nan

melt_y=toplot
        
# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
#pos = ax1.imshow(melt_y, vmin=-2, vmax=+2)
vmin=-2
vmax=+2
cmap = cmap = mpl.cm.PiYG
bounds = np.linspace(vmin,vmax,17)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot, norm=norm, cmap='PiYG')

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title(f'melt year difference [yr]\n' + r'$meltYear(\rho=910kgm^{-3})-meltYear(\rho=850kgm^{-3})$')
ax1.set_aspect(1090.5/1273.5)

plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')

#fig.colorbar(pos, ax=ax1)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1)#norm=norm, cmap=cmap),ax=ax1)
plt.savefig(f'{out_path}/40x40/plots/fig310b.png')
plt.show()

#====================
#i_ens=""
#i_ssp=":"
#dum = np.zeros(n * n)
#dum[g_ind] = ensemble_mean910 - ensemble_mean850
#dum_resh = np.reshape(dum, [n, n])

#melt_y=dum_resh
        
# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
#fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
#pos = ax1.imshow(melt_y, vmin=-2, vmax=2)

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

#ax1.set_title(f'melt year difference (yr)' + r'$meltYear(\rho=910kgm^{-3})-meltYear(\rho=850kgm^{-3})$')
#ax1.set_aspect(1090.5/1273.5)
#plt.ylabel('latitude from origin latitude [m]')
#plt.xlabel('longitude from origin longitude [m]')
#
#fig.colorbar(pos, ax=ax1)
#plt.show()

#Figure 3.11: Melt year for each point independently, with ρice = 850kgm3, H=5m, no
#surface distinction. Ensemble mean for all GCMs and ssp scenarios
ensemble_mean790=ensemble_sum790/15
ensemble_mean850=ensemble_sum850/15
ensemble_mean910=ensemble_sum910/15
ensemble_mean_nosd=ensemble_sum_nosd/15



dum = np.zeros(n * n)
dum[g_ind] = ensemble_mean_nosd - ensemble_mean850
toplot = np.reshape(dum, [n, n])

toplot[tmp_limits] = np.nan

toplot[points_outliers_i] = np.nan

melt_y=toplot
        
# plot
# https://matplotlib.org/stable/gallery/color/colorbar_basics.html
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
# plot just the positive data and save the
# color "mappable" object returned by ax1.imshow
#pos = ax1.imshow(melt_y, vmin=-2, vmax=+2)
vmin=0
vmax=20
cmap = mpl.cm.viridis
bounds = np.linspace(vmin,vmax,17)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

pos = ax1.imshow(toplot, norm=norm, cmap=cmap)

# add the colorbar using the figure's method,
# telling which mappable we're talking about and
# which axes object it should be near

plt.xticks(default_x_ticks, x, rotation=-45)
plt.yticks(default_y_ticks, y)
ax1.set_title(f'melt year [yr]\n' + r'$meltYear(surf.dis=no)-meltYear(surf.dis=yes)$')
ax1.set_aspect(1090.5/1273.5)

plt.ylabel('latitude from origin latitude [m]')
plt.xlabel('longitude from origin longitude [m]')

#fig.colorbar(pos, ax=ax1)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1)#norm=norm, cmap=cmap),ax=ax1)
plt.savefig(f'{out_path}/40x40/plots/fig311.png')
plt.show()

###############################################################################
#Figure 3.12: Glacier area evolution during the 2020-2060 period. ρice = 850kgm3,
#H=5m. Ensemble prediction for all GCMs and ssp scenarios. Colours represent the 3 ssp
#...include legend!...

# area #
fig, ax1 = plt.subplots(figsize=(7, 7), dpi=200,ncols=1)
for i_ens in [0,1,2,3,4]:
    for i_ssp in [0,1,2]:
        if i_ssp==0 :
            ltype='0.05'
        label={i_ens}
        plt.plot(np.array(obj.index[0:412]+2020.68),np.array(proj_pd.iloc[1][i_ens][i_ssp][1][0:412])/proj_pd.iloc[1][i_ens][i_ssp][1][0],
                 ltype,
                 linewidth=0.3
                 )
plt.plot(np.array(obj.index[0:412]+2020.68),np.array(proj_pd.iloc[1][i_ens][i_ssp][1][0:412])/proj_pd.iloc[1][i_ens][i_ssp][1][0],
         ltype,
         linewidth=0.4,
         label=f'projection members'
         )
a00=np.asarray((proj_pd.iloc[1][0][0][1]))/proj_pd.iloc[1][0][2][1][0]
a11=np.asarray((proj_pd.iloc[1][0][0][1]))/proj_pd.iloc[1][0][2][1][0]
a22=np.asarray((proj_pd.iloc[1][0][0][1]))/proj_pd.iloc[1][0][2][1][0]
a03=np.asarray((proj_pd.iloc[1][1][0][1]))/proj_pd.iloc[1][0][2][1][0]
a14=np.asarray((proj_pd.iloc[1][1][0][1]))/proj_pd.iloc[1][0][2][1][0]
a20=np.asarray((proj_pd.iloc[1][1][1][1]))/proj_pd.iloc[1][0][2][1][0]
a01=np.asarray((proj_pd.iloc[1][2][1][1]))/proj_pd.iloc[1][0][2][1][0]
a12=np.asarray((proj_pd.iloc[1][2][1][1]))/proj_pd.iloc[1][0][2][1][0]
a23=np.asarray((proj_pd.iloc[1][2][1][1]))/proj_pd.iloc[1][0][2][1][0]
a04=np.asarray((proj_pd.iloc[1][3][1][1]))/proj_pd.iloc[1][0][2][1][0]
a10=np.asarray((proj_pd.iloc[1][3][2][1]))/proj_pd.iloc[1][0][2][1][0]
a21=np.asarray((proj_pd.iloc[1][3][2][1]))/proj_pd.iloc[1][0][2][1][0]
a02=np.asarray((proj_pd.iloc[1][4][2][1]))/proj_pd.iloc[1][0][2][1][0]
a13=np.asarray((proj_pd.iloc[1][4][2][1]))/proj_pd.iloc[1][0][2][1][0]
a24=np.asarray((proj_pd.iloc[1][4][2][1]))/proj_pd.iloc[1][0][2][1][0]

#mean
aa=(a00+a01+a02+a03+a04+
    a10+a11+a12+a13+a14+
    a20+a21+a22+a23+a24)/15
plt.plot(np.array(obj.index[0:412]+2020.68), 
         aa[0:412],
         'red',
         linewidth=1,
         label=f'ensemble mean'
         )
plt.hlines(0.1, xmin=2020.68,xmax=2055, label='10% threshold')
plt.hlines(0.041, xmin=2020.68,xmax=2055,linestyles='--', label='2 ha threshold')

ax1.set_aspect(25)
plt.title('Relative area change')
plt.ylabel('Area relative to Oct 2020 [-]')
plt.xlabel('calendar year')
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.savefig(f'{out_path}/40x40/plots/fig312.png')
plt.show()


#Figure 3.13: Glacier volume evolution during the 2020-2060 period. ρice = 850kgm3,
#H=5m. Ensemble prediction for all GCMs and ssp scenarios. Colours represent the 3 ssp
#...include legend!..

# Volume #
fig, ax1 = plt.subplots(figsize=(7, 7),dpi=200, ncols=1)
for i_ens in [0,1,2,3,4]:
    for i_ssp in [0,1,2]:
        if i_ssp==0 :
            ltype='0.05'
        label={i_ens}
        plt.plot(np.array(obj.index[0:412]+2020.68),np.array(proj_pd.iloc[1][i_ens][i_ssp][0][0:412])/proj_pd.iloc[1][i_ens][i_ssp][0][0],
                 ltype,
                 linewidth=0.3
                 )
plt.plot(np.array(obj.index[0:412]+2020.68),np.array(proj_pd.iloc[1][i_ens][i_ssp][0][0:412])/proj_pd.iloc[1][i_ens][i_ssp][0][0],
         ltype,
         linewidth=0.4,
         label=f'projection members'
         )
b00=np.asarray((proj_pd.iloc[1][0][0][0]))/proj_pd.iloc[1][0][2][0][0]
b11=np.asarray((proj_pd.iloc[1][0][0][0]))/proj_pd.iloc[1][0][2][0][0]
b22=np.asarray((proj_pd.iloc[1][0][0][0]))/proj_pd.iloc[1][0][2][0][0]
b03=np.asarray((proj_pd.iloc[1][1][0][0]))/proj_pd.iloc[1][0][2][0][0]
b14=np.asarray((proj_pd.iloc[1][1][0][0]))/proj_pd.iloc[1][0][2][0][0]
b20=np.asarray((proj_pd.iloc[1][1][1][0]))/proj_pd.iloc[1][0][2][0][0]
b01=np.asarray((proj_pd.iloc[1][2][1][0]))/proj_pd.iloc[1][0][2][0][0]
b12=np.asarray((proj_pd.iloc[1][2][1][0]))/proj_pd.iloc[1][0][2][0][0]
b23=np.asarray((proj_pd.iloc[1][2][1][0]))/proj_pd.iloc[1][0][2][0][0]
b04=np.asarray((proj_pd.iloc[1][3][1][0]))/proj_pd.iloc[1][0][2][0][0]
b10=np.asarray((proj_pd.iloc[1][3][2][0]))/proj_pd.iloc[1][0][2][0][0]
b21=np.asarray((proj_pd.iloc[1][3][2][0]))/proj_pd.iloc[1][0][2][0][0]
b02=np.asarray((proj_pd.iloc[1][4][2][0]))/proj_pd.iloc[1][0][2][0][0]
b13=np.asarray((proj_pd.iloc[1][4][2][0]))/proj_pd.iloc[1][0][2][0][0]
b24=np.asarray((proj_pd.iloc[1][4][2][0]))/proj_pd.iloc[1][0][2][0][0]

#mean
bb=(b00+b01+b02+b03+b04+
    b10+b11+b12+b13+b14+
    b20+b21+b22+b23+b24)/15
plt.plot(np.array(obj.index[0:412]+2020.68), 
         bb[0:412],
         'red',
         linewidth=1,
         label=f'ensemble mean'
         )
plt.plot(np.array(obj.index[0:412]+2020.68), 
         aa[0:412],
         'r--',
         linewidth=1,
         label=f'ensemble mean (area)'
         )

plt.hlines(0.1, xmin=2020.68,xmax=2055, label=f'10% threshold')
#plt.hlines(0.08, colors='blue', xmin=2020.68,xmax=2055,linestyles='--', label=f'10% threshold (area)')

ax1.set_aspect(22)
plt.title('Relative volume change')
plt.ylabel('Volume relative to Oct 2020 [-]')
plt.xlabel('calendar year')
plt.legend(loc='best', bbox_to_anchor=(1, 0.5))
plt.savefig(f'{out_path}/40x40/plots/fig313.png')
plt.show()








