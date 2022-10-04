#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piont mass balance deconstruction. First I build this and then I will create a class in the 
main script. Start with the code from 
"../docs/aneto_test/use_your_inventory.ipynb"
"""
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

cgidf = gpd.read_file(utils.get_demo_file('cgi2.shp'))
cgidf

cgidf_a = gpd.read_file('/home/francesc/data/aneto_glacier/Contornos/Aneto2011.shp')

rgidf_simple_a = utils.cook_rgidf(cgidf_a, ids=[int('3208')], o1_region='11', o2_region='02', bgndate='2011') #id_suffix aneto glacier
rgidf_simple_a

from oggm import cfg, workflow

cfg.initialize()
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['use_rgi_area'] = False
cfg.PARAMS['use_intersects'] = False
cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-sfc-type', reset=True)
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

cfg.initialize() #logging_level='WARNING'
cfg.PARAMS['use_multiprocessing'] = False
cfg.PARAMS['continue_on_error'] = False
cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-sfc-type')#, reset=True)
cfg.PARAMS['use_intersects'] = False

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

from MBsandbox.mbmod_daily_oneflowline import process_w5e5_data
from projections_bayescalibration import process_isimip_data


# 2011-2019
#process_w5e5_data(gdir, temporal_resol=temporal_resol,
#                  climate_type=baseline_climate) #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.

# 2020 --> 20xx
process_isimip_data(gdir, temporal_resol=temporal_resol, climate_historical_filesuffix='_daily_WFDE5_CRU') #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.

# MB model options: 
# we just use the most complicated ones: 
mb_type = 'mb_real_daily'  # daily climate resolution 
grad_type = 'var_an_cycle'  # a variable temperature lapse rate (cte between the years but changing spatially and between the months)
melt_f_update = 'monthly'  # how often the melt factor is updated to distinguish between different snow / firn ages
melt_f_change = 'neg_exp'  # the way how the melt factor changes over time 
tau_e_fold_yr = 1  # how fast the melt factor of snow converges to the ice melt factor 
kwargs_for_TIModel_Sfc_Type = {'melt_f_update':melt_f_update, 'melt_f_change': melt_f_change, 'tau_e_fold_yr':tau_e_fold_yr}

melt_f = 200
pf = 1.5

###############################################################################
# here the fun starts: this will be converted into TIModel_Sfc_Type_point class.
hbins = np.NaN # Default is np.NaN
interpolation_optim = False # Default is False

tau_e_fold_yr = tau_e_fold_yr
if melt_f_change == 'neg_exp':
    assert tau_e_fold_yr > 0, "tau_e_fold_yr has to be above zero!"
melt_f_change = melt_f_change
assert melt_f_change in ['linear', 'neg_exp'], "melt_f_change has to be either 'linear' or 'neg_exp'"
# ratio of snow melt_f to ice melt_f
melt_f_ratio_snow_to_ice = 0.5 # Default is 0.5
melt_f_update = melt_f_update
spinup_yrs = 6 # Default is 6 years
# amount of bucket and bucket naming depends on melt_f_update:
if melt_f_update == 'monthly':
    # for each month over 6 years a bucket-> in total 72!
    buckets = np.arange(0, 12 * 6, 1).tolist()
else:
    raise Exception('melt_f_update has to be monthly in our case!') #InvalidParamsError('melt_f_update can either be annual or monthly!')
# first bucket: if annual it is 'snow', if monthly update it is 0
first_snow_bucket = buckets[0]

columns = buckets + ['delta_kg/m2']
# TODO: maybe also include snow_delta_kg/m2, firn_yr_1_delta_kg/m2...
#  (only important if I add refreezing or a densification scheme)
# comment: I don't need an ice bucket because this is assumed to be "infinite"
# (instead just have a 'delta_kg/m2' bucket)

# save the inversion height to later check if the same height is applied!!!
inv_heights = surface_h = 3200 # Dummy value for now
check_availability = True # Default is True

# container template (has to be updatable -> pd_mb_template is property/setter thing)
# use the distance_along_flowline as index
_pd_mb_template = pd.DataFrame(0, index=np.arange(0,1),
                                    columns=[]) # exectime-> directly addind columns here should be faster
_pd_mb_template.index.name = 'distance_along_flowline'

# bucket template:
# make a different template for the buckets, because we want directly the right columns inside
_pd_mb_template_bucket = pd.DataFrame(0, index=[3200.4],
                                    columns=columns)  # exectime-> directly addind columns here should be faster
_pd_mb_template_bucket.index.name = 'altitude (masl)'

# storage containers for monthly and annual mb
# columns are the months or years respectively
# IMPORTANT: those need other columns
pd_mb_monthly = _pd_mb_template.copy()
pd_mb_annual = _pd_mb_template.copy()

# bucket containers with buckets as columns
pd_bucket = _pd_mb_template_bucket.copy()
# exectime comment:  this line is quite expensive -> actually this line can be removed as
# _pd_mb_template was defined more clever !!!
# pd_bucket[columns] = 0
# I don't need a total_kg/m2 because I don't know it anyway
# as we do this before inversion!


# storage container for buckets (in kg/m2)
# columns are the buckets
# (6*12 or 6 buckets depending if melt_f_update is monthly or annual)
pd_bucket = pd_bucket

mb_geodetic = np.array(-800) # from zaragoza people 2011-2020

years = np.arange(2011-6, 2040, 1)
ys=2011-6
#ye=2020
ye=2040
#mb_specific = mb_mod_monthly_0_5_m.get_specific_mb(heights=h,
#                                        widths=w,
#                                        year=years
#                                        )#.mean()

# as self.prcp is produced by changing prcp_fac
# we need to update the prcp via prcp_fac by a property (see code after __init__)
# (reason why there is here only a self._prcp_fac and no self.prcp_fac)
#  to allow prcp_fac to be changed after instantiation
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

rho=800

input_filesuffix = '_daily_{}'.format(baseline_climate)
filename = 'climate_historical'

fpath = gdir.get_filepath(filename, filesuffix=input_filesuffix)


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
### not sure
#    if type(self) == TIModel_Sfc_Type:
#        if hydro_month_start != 1 or cfg.PARAMS[f'hydro_month_{self.hemisphere}'] != 1:
#            raise InvalidWorkflowError('TIModel_Sfc_Type works only with calendar years, set '
#                                       'cfg.PARAMS["hydro_month_nh"] (or sh)'
#                                       ' to 1 and process the climate data again')

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
    ref_hgt = xr_nc.ref_hgt # i guess refference hgt is the default altitude of the climate data
    # if climate dataset has been corrected once again
    # or non corrected reference height!
    #try:
    #    self.uncorrected_ref_hgt = xr_nc.uncorrected_ref_hgt
    #except:
    #    self.uncorrected_ref_hgt = xr_nc.ref_hgt

    ys = years[0] if ys is None else ys
    ye = years[-1] if ye is None else ye


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
        temp2d = np.atleast_2d(itemp).repeat(npix, 0) + grad_temp
        # temp_for_melt is computed separately depending on mb_type
        # todo: exectime the line code below is quite expensive in TIModel
        temp2dformelt = _get_tempformelt(temp2d, pok)

        # Compute solid precipitation from total precipitation
        prcp = np.atleast_2d(iprcp).repeat(npix, 0)
        fac = 1 - (temp2d - t_solid) / (t_liq - t_solid)
        # line code below also quite expensive!
        prcpsol = prcp * clip_array(fac, 0, 1)
        return temp2d, temp2dformelt, prcp, prcpsol

####
#get lapse rate and climate data, if will be used in _get_climate:




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
def get_monthly_mb(heights, year=None, add_climate=False,
                   **kwargs):
    """ computes annual mass balance in m of ice per second!
    Attention year is here in hydro float year
    year has to be given as float hydro year from what the month is taken,
    hence year 2000 -> y=2000, m = 1, & year = 2000.09, y=2000, m=2 ...
    which corresponds to the real year 1999 and months October or November
    if hydro year starts in October
    """
    # todo: can actually remove **kwargs???
    # comment: get_monthly_mb and get_annual_mb are only different
    #  to OGGM default for mb_real_daily

    if mb_type == 'mb_real_daily':
        # get 2D values, dependencies on height and time (days)
        out = _get_2d_monthly_climate(heights, year)
        t, temp2dformelt, prcp, prcpsol = out
        # (days per month)
        # dom = 365.25/12  # len(prcpsol.T)
        fact = 12/365.25
        # attention, I should not use the days of years as the melt_f is
        # per month ~mean days of that year 12/daysofyear
        # to have the same unit of melt_f, which is
        # the monthly temperature sensitivity (kg /m² /mth /K),
        mb_daily = prcpsol - melt_f * temp2dformelt * fact

        mb_month = np.sum(mb_daily, axis=1)
        # more correct than using a mean value for days in a month
        #warnings.warn('there might be a problem with SEC_IN_MONTH'
        #              'as February changes amount of days inbetween the years'
        #              ' see test_monthly_glacier_massbalance()')

    # residual is in mm w.e per year, so SEC_IN_MONTH ... but mb_month
    # should be per month!
    mb_month -= residual * SEC_IN_MONTH / SEC_IN_YEAR
    # this is for mb_pseudo_daily otherwise it gives the wrong shape
    mb_month = mb_month.flatten()
    # if add_climate: #default is False
        # if self.mb_type == 'mb_real_daily':
        #     # for run_with_hydro want to get monthly output (sum of daily),
        #     # if we want daily output in run_with_hydro need to directly use get_daily_mb()
        #     prcp = prcp.sum(axis=1)
        #     prcpsol = prcpsol.sum(axis=1)
        #     t = t.mean(axis=1)
        #     temp2dformelt = temp2dformelt.sum(axis=1)
        # if self.mb_type == 'mb_pseudo_daily':
        #     temp2dformelt = temp2dformelt.flatten()
        # return (mb_month / SEC_IN_MONTH / self.rho, t, temp2dformelt,
        #         prcp, prcpsol)
    # instead of SEC_IN_MONTH, use instead len(prcpsol.T)==daysinmonth
    return mb_month / SEC_IN_MONTH / rho

#######
def _get_2d_annual_climate(heights, year):
    return _get_climate(heights, 'annual', year=year)
