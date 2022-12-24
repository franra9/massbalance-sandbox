#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piont mass balance deconstruction. First I build this and then I will create a class in the 
main script. Start with the code from 
"../docs/aneto_test/use_your_inventory.ipynb"
"""
# import parameters
from params import *
import params as params

#  import stuff
import numpy as np
import pandas as pd
import xarray as xr
import warnings
import logging

# Let's get the sample CGI2 glacier inventory and see what it looks like
from oggm import utils
from oggm.utils import (floatyear_to_date, clip_array, clip_min)
import geopandas as gpd


cgidf_a = gpd.read_file('/home/francesc/data/aneto_glacier/Contornos/Aneto2011.shp')

rgidf_simple_a = utils.cook_rgidf(cgidf_a, ids=[int('3208')], o1_region='11', o2_region='02', bgndate='2011') #id_suffix aneto glacier
#rgidf_simple_a

from oggm import cfg, workflow

cfg.initialize()
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['use_rgi_area'] = False
cfg.PARAMS['use_intersects'] = False
cfg.PATHS['working_dir'] = '/home/francesc/data/aneto_glacier/climate/'#utils.gettempdir(dirname='OGGM-sfc-type')#, reset=True) # TODO: mydir
gdirs = workflow.init_glacier_directories(rgidf_simple_a)

# The tasks below require downloading new data - we comment them for the tutorial, but it should work for you!
# workflow.gis_prepro_tasks(gdirs)
# workflow.download_ref_tstars('https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5')
# workflow.climate_tasks(gdirs)
# workflow.inversion_tasks(gdirs)

utils.compile_glacier_statistics(gdirs)

rgidf_simple_a.to_crs('EPSG:25831').area #in m²
rgidf_simple_a['Area'] = rgidf_simple_a.to_crs('EPSG:25831').area * 1e-6 #in km²
rgidf_simple_a
#gdirs

log = logging.getLogger(__name__)

#cfg.initialize() #logging_level='WARNING'


# use Huss flowlines # I don't use this
base_url = ('https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/'
            'L1-L2_files/elev_bands')

# FRA: WHY GEODETIC MB STARTS ON JAN AND NOT SEP?
cfg.PARAMS['hydro_month_nh']=1  # to calibrate to the geodetic estimates we need calendar years !!!

#df = ['RGI60-11.03208']  # here we select Aneto glacier!

gdirs = workflow.init_glacier_directories(rgidf_simple_a,prepro_base_url=base_url)
gdir = gdirs[0]
gdir

temporal_resol = 'daily'
baseline_climate = 'W5E5'
#baseline_climate = 'ISIMIP3B'
#baseline_climate = 'WFDE5_CRU'
#baseline_climate = 'WFDE5_CRU'


#from MBsandbox.mbmod_daily_oneflowline import process_w5e5_data
#from projections_bayescalibration import process_isimip_data
#from MBsandbox.wip.projections_bayescalibration import process_isimip_data

# 2011-2019
#if cal == 'w5e5': 
#process_w5e5_data(gdir, temporal_resol=temporal_resol,
#                  climate_type=baseline_climate) #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.

# 2020 --> 20xx
#if cal == 'isimip3b':
#    process_w5e5_data(gdir, temporal_resol=temporal_resol,
#                      climate_type=baseline_climate) #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.
# TODO: do loop here
#    ssp = params.ssp

#for ssp in ssps:
#    for ensamble_name in ensamble_names[:-1]:
#        process_isimip_data(gdir,
#                        temporal_resol=temporal_resol,
#                        climate_historical_filesuffix='_daily_W5E5',
#                        #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.
#                        #climate_historical_filesuffix = '',
#                        ensemble = ensamble_name,#'mri-esm2-0_r1i1p1f1',
#                        # from temperature tie series the "median" ensemble
#                        ssp = ssp, flat = True,
#                        cluster = False,
#                        year_range = ('2005', '2099'),
#                        correct = False)

# 245 file does not exist!
#process_isimip_data(gdir,
#                    temporal_resol=temporal_resol,
#                    climate_historical_filesuffix='_daily_WFDE5_CRU',
#                    #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.
#                    #climate_historical_filesuffix = '',
#                    ensemble = 'mri-esm2-0_r1i1p1f1',
#                    # from temperature tie series the "median" ensemble
#                    ssp = 'ssp245', flat = True,
#                    cluster = False,
#                    year_range = ('2005', '2050'),
#                    correct = False)

#process_isimip_data(gdir,
#                    temporal_resol=temporal_resol,
#                    climate_historical_filesuffix='_daily_WFDE5_CRU',
#                    #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.
#                    #climate_historical_filesuffix = '',
#                    ensemble = 'mri-esm2-0_r1i1p1f1',
#                    # from temperature tie series the "median" ensemble
#                    ssp = 'ssp370', flat = True,
#                    cluster = False,
#                    year_range = ('2005', '2050'),
#                    correct = False)

#process_isimip_data(gdir,
#                    temporal_resol=temporal_resol,
#                    climate_historical_filesuffix='_daily_WFDE5_CRU',
#                    #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.
#                    #climate_historical_filesuffix = '',
#                    ensemble = 'mri-esm2-0_r1i1p1f1',
#                    # from temperature tie series the "median" ensemble
#                    ssp = 'ssp585', flat = True,
#                    cluster = False,
#                    year_range = ('2005', '2050'),
#                    correct = False)

# MB model options: 
# we just use the most complicated ones: 
mb_type = 'mb_real_daily'  # daily climate resolution 
grad_type = 'var_an_cycle'  # a variable temperature lapse rate (cte between the years but changing spatially and between the months)
melt_f_update = 'monthly'  # how often the melt factor is updated to distinguish between different snow / firn ages
melt_f_change = 'neg_exp'  # the way how the melt factor changes over time 
tau_e_fold_yr = 1  # how fast the melt factor of snow converges to the ice melt factor 
kwargs_for_TIModel_Sfc_Type = {'melt_f_update':melt_f_update, 'melt_f_change': melt_f_change, 'tau_e_fold_yr':tau_e_fold_yr}

###############################################################################
#  prescribe the prcp_fac as it is instantiated
prcp_fac = 2.89
_prcp_fac = prcp_fac #because it is takes as constant.
# same is true for temp bias
_temp_bias = 0.
residual = 0.
default_grad = -0.0065 #K/m, lapse rate
temp_local_gradient_bounds = [-0.009, -0.003]

# temperature threshold where snow/ice melts
t_melt = 0.
t_solid = 0
t_liq = 2

if cal == 'w5e5':
    #input_filesuffix = '_daily_{}'.format(baseline_climate)
    #filename = 'climate_historical'
    #fpath = gdir.get_filepath(filename, filesuffix=input_filesuffix)
    fpath = '/home/francesc/data/aneto_glacier/climate/per_glacier/RGI60-11/RGI60-11.03/RGI60-11.03208/climate_historical_daily_W5E5.nc'
if cal == 'isimip3b':
    fpath = f"/home/francesc/data/aneto_glacier/climate/per_glacier/RGI60-11/RGI60-11.03/RGI60-11.03208/gcm_data_daily_ISIMIP3b_{ensamble_name}_{ssp}_no_correction.nc"

#xr.open_dataset(fpath) as xr_nc
with xr.open_dataset(fpath) as xr_nc:
    if mb_type == 'mb_real_daily' or mb_type == 'mb_monthly':
        # even if there is temp_std inside the dataset, we won't use
        # it for these mb_types
        temp_std = np.NaN

    # goal is to get self.years/self.months in hydro_years
    # (only important for TIModel if we do not use calendar years,
    # for TIModel_Sfc_Type, we can only use calendar years anyways!)

    # get the month where the hydrological month starts
    # as chosen from the gdir climate file
    hydro_month_start = int(xr_nc.time[0].dt.month.values)
    # if we actually use TIModel_Sfc_Type -> hydro_month has to be 1

    if mb_type == 'mb_real_daily':
        # use pandas to convert month/year to hydro_years
        # this has to be done differently than above because not
        # every month, year has the same amount of days
        pd_test = pd.DataFrame(xr_nc.time.to_series().dt.year.values,
                               columns=['year'])
        pd_test.index = xr_nc.time.to_series().values
        pd_test['month'] = xr_nc.time.to_series().dt.month.values
        pd_test['hydro_year'] = np.NaN

        if hydro_month_start == 1:
            # hydro_year corresponds to normal year
            pd_test.loc[pd_test.index.month >= hydro_month_start,
                        'hydro_year'] = pd_test['year']
        else:
            pd_test.loc[pd_test.index.month < hydro_month_start,
                        'hydro_year'] = pd_test['year']
            # otherwise, those days with a month>=hydro_month_start
            # belong to the next hydro_year
            pd_test.loc[pd_test.index.month >= hydro_month_start,
                        'hydro_year'] = pd_test['year']+1
        # month_hydro is 1 if it is hydro_month_start
        month_hydro = pd_test['month'].values+(12-hydro_month_start+1)
        month_hydro[month_hydro > 12] += -12
        pd_test['hydro_month'] = month_hydro
        pd_test = pd_test.astype('int')
        years = pd_test['hydro_year'].values
        ny = years[-1] - years[0]+1
        months = pd_test['hydro_month'].values
    # Read timeseries and correct it
    temp = xr_nc['temp'].values.astype(np.float64) + _temp_bias
    # this is prcp computed by instantiation
    # this changes if prcp_fac is updated (see @property)
    prcp = xr_nc['prcp'].values.astype(np.float64) * _prcp_fac
    print (years)
    # lapse rate (temperature gradient)
    if grad_type == 'var' or grad_type == 'var_an_cycle':
        try:
            # need this to ensure that gradients are not fill-values
            xr_nc['gradient'] = xr_nc['gradient'].where(xr_nc['gradient'] < 1e12)
            grad = xr_nc['gradient'].values.astype(np.float64)
            # Security for stuff that can happen with local gradients
            g_minmax = temp_local_gradient_bounds

            # if gradient is not a number, or positive/negative
            # infinity, use the default gradient
            grad = np.where(~np.isfinite(grad), default_grad, grad)

            # if outside boundaries of default -0.009 and above
            # -0.003 -> use the boundaries instead
            grad = clip_array(grad, g_minmax[0], g_minmax[1])

            if grad_type == 'var_an_cycle':
                # if we want constant lapse rates over the years
                # that change over the annual cycle, but not over time
                if mb_type == 'mb_real_daily':
                    grad_gb = xr_nc['gradient'].groupby('time.month')
                    grad = grad_gb.mean().values
                    g_minmax = temp_local_gradient_bounds

                    # if gradient is not a number, or positive/negative
                    # infinity, use the default gradient
                    grad = np.where(~np.isfinite(grad), default_grad,
                                    grad)
                    assert np.all(grad < 1e12)
                    # if outside boundaries of default -0.009 and above
                    # -0.003 -> use the boundaries instead
                    grad = clip_array(grad, g_minmax[0], g_minmax[1])

                    stack_grad = grad.reshape(-1, 12)
                    grad = np.tile(stack_grad.mean(axis=0), ny)
                    reps_day1 = xr_nc.time[xr_nc.time.dt.day == 1]
                    reps = reps_day1.dt.daysinmonth
                    grad = np.repeat(grad, reps)

                else:
                    stack_grad = grad.reshape(-1, 12)
                    grad = np.tile(stack_grad.mean(axis=0), ny)
        except KeyError:
            text = ('there is no gradient available in chosen climate'
                    'file, try instead e.g. ERA5_daily or ERA5dr e.g.'
                    'oggm.shop.ecmwf.process_ecmwf_data'
                    '(gd, dataset="ERA5dr")')

            raise InvalidParamsError(text)

    elif grad_type == 'cte':
        # if grad_type is chosen cte, we use the default_grad!
        grad = prcp * 0 + default_grad
    else:
        raise InvalidParamsError('grad_type can be either cte,'
                                 'var or var_an_cycle')
    grad = grad
    ref_hgt = xr_nc.ref_hgt 
    # if climate dataset has been corrected once again
    # or non corrected reference height!
    #try:
    #    self.uncorrected_ref_hgt = xr_nc.uncorrected_ref_hgt
    #except:
    #    self.uncorrected_ref_hgt = xr_nc.ref_hgt

#    ys = spin_years[0] if ys is None else ys
#    ye = years[-1] if ye is None else ye


####
#
#######
def _get_tempformelt(temp, pok):
    """ Helper function to compute tempformelt to avoid code duplication
    in get_monthly_climate() and _get2d_annual_climate()
    comment: I can't use  _get_tempformelt outside the class, but sometimes this could be useful.
    If using this again outside of this class, need to remove the "self",
    such as for 'mb_climate_on_height' in climate.py, that has no self....
    (would need to change temp, t_melt ,temp_std, mb_type, N, loop)
    Parameters
    -------
        temp: ndarray
            temperature time series
        pok: ndarray
            indices of time series
    Returns
    -------
    (tempformelt)
    """

    tempformelt_without_std = temp - t_melt

    # computations change only if 'mb_pseudo_daily' as mb_type!
    if mb_type == 'mb_monthly' or mb_type == 'mb_real_daily':
        tempformelt = tempformelt_without_std

    else:
        raise InvalidParamsError('mb_type can only be "mb_monthly,\
                                 mb_pseudo_daily or mb_real_daily" ')
    #  replace all values below zero to zero
    # todo: exectime this is also quite expensive
    clip_min(tempformelt, 0, out=tempformelt)

    return tempformelt



####
# define here _get_climate?
def _get_climate(heights, climate_type, year=None):
    """ Climate information at given heights.
    Note that prcp is corrected with the precipitation factor and that
    all other model biases (temp and prcp) are applied.
    same as in OGGM default except that tempformelt is computed by
    self._get_tempformelt
    Parameters
    -------
    heights : np.array or list
        heights along flowline
    climate_type : str
        either 'monthly' or 'annual', if annual floor of year is used,
        if monthly float year is converted into month and year
    year : float
        float hydro year from what both, year and month, are taken if climate_type is monthly.
        hence year 2000 -> y=2000, m = 1, & year = 2000.09, y=2000, m=2 ...
        which corresponds to the real year 1999 and months October or November
        if hydro year starts in October
    Returns
    -------
    (temp, tempformelt, prcp, prcpsol)
    """

    y, m = floatyear_to_date(year)
    ys=2005 #temp
    ye=2099
    
    if y < ys or y > ye:
        raise ValueError('year {} out of the valid time bounds: '
                         '[{}, {}]'.format(y, ys, ye))

    if mb_type == 'mb_real_daily' or climate_type == 'annual':
        #print(years)
        pok = np.where((years == y) & (months == m))[0]
        #print(pok)
        if len(pok) < 28:
            warnings.warn('something goes wrong with amount of entries\
                          per month for mb_real_daily')
    # Read time series
    # (already temperature bias and precipitation factor corrected!)
    itemp = temp[pok]
    prcp = xr_nc['prcp'].values.astype(np.float64) * _prcp_fac # problems with the scope if I dont define it here...
    iprcp = prcp[pok]
    igrad = grad[pok]

    # For each height pixel:
    # Compute temp and tempformelt (temperature above melting threshold)
    heights = np.asarray(heights)
    #npix = len(heights)
    npix = 1
    if mb_type == 'mb_real_daily' or climate_type == 'annual':
        grad_temp = np.atleast_2d(igrad).repeat(npix, 0)
        # todo: exectime the line code below is quite expensive in TIModel
        grad_temp *= (heights.repeat(len(pok)).reshape(grad_temp.shape) -
                      ref_hgt)
        #print(f'ref_height: {ref_hgt}, cal:{cal}, grad_temp:{grad_temp}')
        temp2d = np.atleast_2d(itemp).repeat(npix, 0) + grad_temp
        # temp_for_melt is computed separately depending on mb_type
        # todo: exectime the line code below is quite expensive in TIModel
        temp2dformelt = _get_tempformelt(temp2d, pok)

        # Compute solid precipitation from total precipitation
        prcp = np.atleast_2d(iprcp).repeat(npix, 0)
        fac = 1 - (temp2d - t_solid) / (t_liq - t_solid)
        # line code below also quite expensive!
        prcpsol = prcp * clip_array(fac, 0, 1)
        #print(f'Grad_temp {grad_temp.mean()}, data={cal}') 
        #print(f'temp:{temp2d.mean()}: data={cal}')
        #print(f'temp4melt:{temp2dformelt.mean()}: data={cal}')
        #print(f'ref_height {ref_hgt}') 
        #print(f'temp:{temp2d.mean()}, ref_height: {ref_hgt}, cal:{cal}, grad_temp:{grad_temp.mean()}')
        #print(f'cal:{cal}, temp2d:{temp2d.mean()}-itemp:{itemp.mean()}, grad_temp:{grad_temp.mean()}')
        return temp2d, temp2dformelt, prcp, prcpsol

########
# def _get_2d_monthly_climate(self, heights, year=None):
def _get_2d_monthly_climate(heights, year=None): # this is the same as _get climate!
    # only works with real daily climate data!
    # (is used by get_monthly_mb of TIModel)
    # comment: as it is not used in TIModel_Sfc_type it could also directly go inside of TIModel ?!
    if mb_type == 'mb_real_daily':
        return _get_climate(heights, 'monthly', year=year)
    else:
        raise InvalidParamsError('_get_2d_monthly_climate works only\
                                 with mb_real_daily as mb_type!!!')

########
# from # class TIModel(TIModel_Parent):
#def get_monthly_mb(heights, year=None, add_climate=False,
#                   **kwargs):
#    """ computes annual mass balance in m of ice per second!
#    Attention year is here in hydro float year
#    year has to be given as float hydro year from what the month is taken,
#    hence year 2000 -> y=2000, m = 1, & year = 2000.09, y=2000, m=2 ...
#    which corresponds to the real year 1999 and months October or November
#    if hydro year starts in October
#    """
#    # todo: can actually remove **kwargs???
#    # comment: get_monthly_mb and get_annual_mb are only different
#    #  to OGGM default for mb_real_daily
#
#    if mb_type == 'mb_real_daily':
#        # get 2D values, dependencies on height and time (days)
#        out = _get_2d_monthly_climate(heights, year)
#        t, temp2dformelt, prcp, prcpsol = out
#        # (days per month)
#        # dom = 365.25/12  # len(prcpsol.T)
#        fact = 12/365.25
#        # attention, I should not use the days of years as the melt_f is
#        # per month ~mean days of that year 12/daysofyear
#        # to have the same unit of melt_f, which is
#        # the monthly temperature sensitivity (kg /m² /mth /K),
#        mb_daily = prcpsol - melt_f * temp2dformelt * fact
#
#        mb_month = np.sum(mb_daily, axis=1)
#        # more correct than using a mean value for days in a month
#        #warnings.warn('there might be a problem with SEC_IN_MONTH'
#        #              'as February changes amount of days inbetween the years'
#        #              ' see test_monthly_glacier_massbalance()')
#
#    # residual is in mm w.e per year, so SEC_IN_MONTH ... but mb_month
#    # should be per month!
#    mb_month -= residual * SEC_IN_MONTH / SEC_IN_YEAR
#    # this is for mb_pseudo_daily otherwise it gives the wrong shape
#    mb_month = mb_month.flatten()
#    # if add_climate: #default is False
#        # if self.mb_type == 'mb_real_daily':
#        #     # for run_with_hydro want to get monthly output (sum of daily),
#        #     # if we want daily output in run_with_hydro need to directly use get_daily_mb()
#        #     prcp = prcp.sum(axis=1)
#        #     prcpsol = prcpsol.sum(axis=1)
#        #     t = t.mean(axis=1)
#        #     temp2dformelt = temp2dformelt.sum(axis=1)
#        # if self.mb_type == 'mb_pseudo_daily':
#        #     temp2dformelt = temp2dformelt.flatten()
#        # return (mb_month / SEC_IN_MONTH / self.rho, t, temp2dformelt,
#        #         prcp, prcpsol)
#    # instead of SEC_IN_MONTH, use instead len(prcpsol.T)==daysinmonth
#    return mb_month / SEC_IN_MONTH / rho

#######
def _get_2d_annual_climate(heights, year):
    return _get_climate(heights, 'annual', year=year)
