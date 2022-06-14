#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piont mass balance deconstruction. First I build this and then I will create a class in the 
main script. Start with the code from 
"../docs/aneto_test/use_your_inventory.ipynb"
"""
#  import stuff
import skimage
import numpy as np
import os

import pandas as pd
import xarray as xr
import seaborn as sns
import pickle

import matplotlib.pyplot as plt
import matplotlib

from numpy.testing import assert_allclose
import scipy
import scipy.stats as stats
from IPython.core.pylabtools import figsize
import os

import oggm
from oggm import cfg, utils, workflow, tasks, graphics, entity_task
from oggm.utils import date_to_floatyear
from oggm.shop import gcm_climate
from oggm.core import massbalance, flowline, climate
from oggm.cfg import SEC_IN_YEAR, SEC_IN_MONTH, SEC_IN_DAY

import logging

# Let's get the sample CGI2 glacier inventory and see what it looks like
from oggm import utils
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

from MBsandbox.mbmod_daily_oneflowline import process_w5e5_data

process_w5e5_data(gdir, temporal_resol=temporal_resol,
                  climate_type=baseline_climate) #Processes and writes the WFDE5_CRU & W5E5 daily baseline climate data for a glacier.

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
# here the fun starts: this will be onverted into TIModel_Sfc_Type_point class.
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
if melt_f_update == 'annual':
    buckets = ['snow', 'firn_yr_1', 'firn_yr_2', 'firn_yr_3',
                    'firn_yr_4', 'firn_yr_5']
elif melt_f_update == 'monthly':
    # for each month over 6 years a bucket-> in total 72!
    buckets = np.arange(0, 12 * 6, 1).tolist()
else:
    raise Exception('melt_f_update can either be annual or monthly!') #InvalidParamsError('melt_f_update can either be annual or monthly!')
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

years = np.arange(2011, 2020, 1)
    
#mb_specific = mb_mod_monthly_0_5_m.get_specific_mb(heights=h,
#                                        widths=w,
#                                        year=years
#                                        )#.mean()

rho=800


####
#
####

########
# define here _get_climate?
#
#
#######
def _get_2d_annual_climate(heights, year):
    return _get_climate(heights, 'annual', year=year)
    
#### this is adapted from a method from the class
#   mbmod_daily_oneflowline.TIModel_Sfc_Type
####
def get_annual_mb(heights, year=None, unit='m_of_ice',
                  bucket_output=False, spinup=True,
                  add_climate=False,
                  auto_spinup=True,
                  **kwargs):
    """
    computes annual mass balance in m of ice per second

    Parameters
    ----------
    heights : np.array
        at the moment works only with inversion heights!
    year: int
        integer CALENDAR year! if melt_f_update='monthly', it will loop over each month.
    unit : str
        default is 'm of ice', nothing else implemented at the moment!
        comment: include option of metre of glacier where the different densities
         are taken into account but in this case would need to add further columns in pd_buckets
         like: snow_delta_kg/m2 ... and so on (only important if I add refreezing or a densification scheme)
    bucket_output: bool (default is False)
        if True, returns as second output the pd.Dataframe with details about
        amount of kg/m2 for each bucket and height grid point,
        set it to True to visualize how the buckets change over time or for testing
        (these are the buckets before they got updated for the next year!)
    spinup : bool (default is True)
        if a spinup is applied to fill up sfc type buckets beforehand
        (default are 6 years, check under spinup_yrs)
    add_climate: bool (default is False)
        for run_with_hydro (not yet implemented!)
        todo: implement and test it
    auto_spinup: bool (default is true)
        if True, it automatically computes the spinup years (default is 6) beforehand (however in this case,
        these 6 years at the beginning, although saved in pd_mb_annual, had no spinup...)
        todo: maybe need to add a '_no_spinup' to those that had no spinup?
         or save them in a separate pd.DataFrame?
    **kwargs:
        other stuff passed to get_monthly_mb or to _add_delta_mb_vary_melt_f
        **kwargs necessary to take stuff we don't use (like fls...)

    """

    # when we set spinup_yrs to zero, then there should be no spinup occurring, even if spinup
    # is set to True
    if spinup_yrs == 0:
        spinup = False

    if len(pd_mb_annual.columns) > 0:
        if pd_mb_annual.columns[0] > pd_mb_annual.columns[-1]:
            raise Exception('need to run years in ascending order! Maybe reset first!')
    # dirty check that the given heights are in the right order
    assert heights[0] > heights[-1], "heights should be in descending order!"
    # otherwise pd_buckets does not work ...
    # below would be a more robust test, but this is not compatible to all use cases:
    # np.testing.assert_allclose(heights, inv_heights,
    #                           err_msg='the heights should correspond to the inversion heights',
    #                           )
    # np.testing.assert_allclose(heights, mod_heights)

    # we just convert to integer without checking ...
    year = int(year)
    # comment: there was some reason why I commented that code again
    # if not isinstance(year, int):
    #    raise InvalidParamsError('Year has to be the full year for get_annual_mb,'
    #                             'year needs to be an integer')
    if year < 1979+spinup_yrs:
        # most climatic data starts in 1979, so if we want to
        # get the first 6 years can not use a spinup!!!
        # (or would need to think about sth. else, but not so
        # important right now!!!)
        spinup = False
    if year in pd_mb_annual.columns and check_availability:
        # print('takes existing annual mb')
        # if that year has already been computed, and we did not change any parameter settings
        # just get the annual_mb without redoing all the computations
        mb_annual = pd_mb_annual[year].values
        if bucket_output:
            raise Exception('if you want to output the buckets, you need to do'
                                       'reset_pd_mb_bucket() and rerun')

        if add_climate:
            t, temp2dformelt, prcp, prcpsol = _get_2d_annual_climate(heights,
                                                                          year) #=_get_climate(heights, 'annual', year=year)
            return (mb_annual, t.mean(axis=1), temp2dformelt.sum(axis=1),
                    prcp.sum(axis=1), prcpsol.sum(axis=1))
            # raise NotImplementedError('TODO: add_climate has to be implemented!')
        else:
            return mb_annual
    else:
        # do we need to do the spinup beforehand
        # if any of the default 6 spinup years before was not computed,
        # we need to reset all and do the spinup properly -> (i.e. condi = True)
        if melt_f_update == 'annual':
            # so we really need to check if every year exists!!!
            condis = []
            for bef in np.arange(1, spinup_yrs, 1):
                condis.append(int(year - bef) not in pd_mb_annual.columns)
            condi = np.any(np.array(condis))
        elif melt_f_update == 'monthly':
            try:
                # check if first and last month of each spinup years before exists
                condis = []
                for bef in np.arange(1, spinup_yrs, 1):
                    condis.append(date_to_floatyear(year - bef, 1) not in pd_mb_monthly.columns)
                    condis.append(date_to_floatyear(year - bef, 12) not in pd_mb_monthly.columns)
                condi = np.any(np.array(condis))
            except:
                condi = True

        if condi and spinup and auto_spinup:
            # reset and do the spinup:
            # comment: I think we should always reset when
            # doing the spinup and not having computed year==2000
            # problem: the years before year should not be saved up (in get_annual_mb)!
            # (because they are not computed right!!! (they don't have a spinup)
            reset_pd_mb_bucket()
            for yr in np.arange(year-spinup_yrs, year):
                get_annual_mb(heights, year=yr, unit=unit,
                                   bucket_output=False,
                                   spinup=False, add_climate=False,
                                   auto_spinup=False,
                                   **kwargs)

        if spinup:
            # check if the spinup years had been computed (should be inside of pd_mb_annual)
            # (it should have been computed automatically if auto_spinup=True)
            for bef in np.arange(1, 6, 1):
                if int(year-bef) not in pd_mb_annual.columns.values:
                    raise Exception('need to do get_annual_mb of all spinup years'
                                               '(default is 6) beforehand')
        if melt_f_update == 'annual':
            pd_bucket = _add_delta_mb_vary_melt_f(heights, year=year,
                                                            climate_resol='annual')
            mb_annual = ((pd_bucket['delta_kg/m2'].values
                          - residual) / SEC_IN_YEAR / rho)
            if bucket_output:
                # copy because we want to output the bucket that is not yet updated!!!
                pd_bucket = pd_bucket.copy()
            # update to one year later ... (i.e. put the snow / firn into the next older bucket)
            _update()
            # save the annual mb
            # todo Fabi exectime: --> would need to restructure code to remove pd stuff
            pd_mb_annual[year] = mb_annual
            # this is done already if melt_f_update is monthly (see get_monthly_mb for m == 12)
        elif melt_f_update == 'monthly':
            # will be summed up over each month by getting pd_mb_annual
            # that is set in december (month=12)
            for m in np.arange(1, 13, 1):
                floatyear = date_to_floatyear(year, m)
                out = get_monthly_mb(heights, year=floatyear, unit=unit,
                                          bucket_output=bucket_output, spinup=spinup,
                                          add_climate=add_climate,
                                          auto_spinup=auto_spinup,
                                          **kwargs)
                if bucket_output and m == 12:
                    pd_bucket = out[1]
            # get mb_annual that is produced by
            mb_annual = pd_mb_annual[year].values
        if bucket_output:
            return mb_annual, pd_bucket
        else:
            if add_climate:
                t, temp2dformelt, prcp, prcpsol = _get_2d_annual_climate(heights,
                                                                              year)
                return (mb_annual, t.mean(axis=1), temp2dformelt.sum(axis=1),
                        prcp.sum(axis=1), prcpsol.sum(axis=1))
                #raise NotImplementedError('TODO: add_climate has to be implemented!')
            else:
                return mb_annual

        #todo
        # if add_climate:
        #    return (mb_annual, t.mean(axis=1), tmelt.sum(axis=1),
        #            prcp.sum(axis=1), prcpsol.sum(axis=1))

#### this is adapted from a method from the class
#   mbmod_daily_oneflowline.MultipleFlowlineMassBalance_TIModel 
####
def get_specific_mb( heights=None, widths=None, fls=None,
                    year=None, **kwargs):

    """ computes specific mass-balance for each year in [kg /m2]"""

#    if heights is not None or widths is not None:
#        raise ValueError('`heights` and `widths` kwargs do not work with '
#                         'MultipleFlowlineMassBalance!')

    if fls is None:
        fls = fls

# not sure what does this
#    if len(np.atleast_1d(year)) > 1:
#        out = [get_specific_mb(fls=fls, year=yr, **kwargs) for yr in year] # call a method inside the method???
#        return np.asarray(out)

    mbs = []
    widths = []

    mb = get_annual_mb(fl.surface_h, year=year, fls=fls,
                              fl_id=i, **kwargs)
    mbs = np.append(mbs, mb * SEC_IN_YEAR * rho)
    return np.average(mbs, weights=widths)

##### This is what I want at the end:
#get_specific_mb(heights=h, widths=w, year=years).mean()