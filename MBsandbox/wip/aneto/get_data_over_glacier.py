#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preprocess rasters (fill NA values) + cut outline shape 

@author: franra9 12092022
"""
import os
import rioxarray as rio
import shapely.geometry as shpg
import numpy as np
from scipy import interpolate
import geopandas as gpd
    
def get_raster_data(n):
    """
    Imput
    n : resolution; n x n grid.

    Returns
    -------
    cal_pts : list
        list with the raster data for thickness change, empty melt factor, thickness in 2020
        and subglacial topography(altitude).

    """
    # mother path:
    data_path = '/home/francesc/data/aneto_glacier/'
    
    # data path: subglac_altimetry
    dpth_alti = os.path.join(data_path, 'TopografiaSubglaciar/TopoSubglaciar.tif')            
    
    # data path: thickness change 2011-2020
    dpth_thi1120 = os.path.join(data_path, 'CambiosEspesor/Aneto20112020_int.tif')    
    
    # data path: thickness 2020
    dpth_thi20 = os.path.join(data_path, 'EspesorAneto2020/InterpolacionEspsores2020.tif')  
    
    # data path: outlines 2020
    dpth_outl20 = os.path.join(data_path, 'Contornos/Aneto2020.shp')  
    
    #Read topo, thickness change, gpr 2020 
    alti = rio.open_rasterio(dpth_alti)
    thi1120 = rio.open_rasterio(dpth_thi1120)
    thi20 = rio.open_rasterio(dpth_thi20)

    ### interpolate fo fill missing values
    array = alti[0]
    
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    #mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                              (xx, yy),
                                 method='nearest')
    
    alti[0]=GD1
    
    ### interpolate fo fill missing values
    # thickness change
    array = thi1120[0]
    
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    #mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                              (xx, yy),
                                 method='linear')
    thi1120[0]=GD1
    

    # load glacier outlines
    outlines = gpd.read_file(os.path.join(data_path, dpth_outl20))

    # cut in outline shape
    # select 0-th glacier
    crp0 = outlines.iloc[[0]]
    
    # select 1-th glacier
    crp1 = outlines.iloc[[1]]
    
    # crop
    thi1120_0_cut = thi1120.rio.clip(crp0.buffer(1).apply(shpg.mapping), 
                               outlines.crs, drop=True)
    
    thi1120_1_cut = thi1120.rio.clip(crp1.buffer(1).apply(shpg.mapping), 
                               outlines.crs, drop=True)
    
    # crop
    alti0_cut = alti.rio.clip(crp0.buffer(1).apply(shpg.mapping), 
                               outlines.crs, drop=True)

    alti1_cut = alti.rio.clip(crp1.buffer(1).apply(shpg.mapping),
                               outlines.crs, drop=True)

    # crop
    thi200_cut = thi20.rio.clip(crp0.buffer(1).apply(shpg.mapping), 
                               outlines.crs, drop=True)
    
    thi201_cut = thi20.rio.clip(crp1.buffer(1).apply(shpg.mapping), 
                               outlines.crs, drop=True)

    # regrid to coarser grid (thickness 2000)
    thi1120_0_cut = thi1120_0_cut[0][np.linspace(0,727, thi200_cut.shape[1], dtype=int), np.round(np.linspace(0,849, thi200_cut.shape[2],  dtype=int),0)]
    alti0_cut = alti0_cut[0][np.linspace(0,1090, thi200_cut.shape[1], dtype=int), np.round(np.linspace(0,1275, thi200_cut.shape[2],  dtype=int),0)]

    #thi1120_1_cut = thi1120_1_cut[0][np.linspace(0,138, thi201_cut.shape[1], dtype=int), np.round(np.linspace(0,172, thi201_cut.shape[2],  dtype=int),0)]
    #alti1_cut = alti1_cut[0][np.linspace(0,208, thi201_cut.shape[1], dtype=int), np.round(np.linspace(0,258, thi201_cut.shape[2],  dtype=int),0)]
    
    # select points
    cal_pts = []
    # n  --> comes from find_mf
    
    cal_ind = [np.linspace(0, 362, n, dtype=int), np.linspace(0, 424, n, dtype=int)]
    
    #store as output list
    cal_pts.append(['th1120', 'melt_f', 'th2020', 'alti']) # names
    cal_pts.append(thi1120_0_cut[cal_ind[0],cal_ind[1]])
    cal_pts.append(thi1120_0_cut[cal_ind[0],cal_ind[1]] * 0) # inicialize melt_f pd.df
    cal_pts.append(thi200_cut[0][cal_ind[0],cal_ind[1]])
    cal_pts.append(alti0_cut[cal_ind[0],cal_ind[1]])

    print(f'Values from different imput raster can be found now under the "cal_pts" variable with the following order: {cal_pts[0]}' )
    
    return cal_pts




