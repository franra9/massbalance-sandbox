#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:18:00 2022

@author: francesc
"""
import sys
import logging

print('Number of arguments:'), print(len(sys.argv), 'arguments.')
print('Argument List:'), print(str(sys.argv))

# thing to run many simulations at a time
#ens = int(sys.argv[1])
#ss = int(sys.argv[2])
#n = int(sys.argv[3])
#rho = int(sys.argv[4])

import numpy as np

# calibration period: 2011-2019 (w5e5) # i should change thisvariable. for projection is not useful anymore.

cal = 'w5e5'
#cal = 'isimip3b'

ensamble_names = ['ukesm1-0-ll_r1i1p1f2', 'gfdl-esm4_r1i1p1f1', 'ipsl-cm6a-lr_r1i1p1f1',
       'mpi-esm1-2-hr_r1i1p1f1', 'mri-esm2-0_r1i1p1f1', 'w5e5']

ens=5 #2
ss=1 #1
n=40#40
rho=850

rho = rho/1000
#rho = 0.9 # 900kg/m3

ensamble_name=ensamble_names[ens]
#ensamble_name='gfdl-esm4_r1i1p1f1'

#'ssp126' #'ssp126' #"ssp245" does not exist!
ssps = ['ssp126', 'ssp370', 'ssp585']
ssp = ssps[ss]
#ssp = 'ssp370' 
#ssp = 'ssp585'
# paralelization:
ncores = 5

if cal == 'w5e5': # do some assertion here
    ssp = ''
    y_alfa = 2011
    y_omega = 2019
    ensamble_name = str(cal)
    
if cal == 'isimip3b': # do some assertion here
    y_alfa = 2011
    y_omega = 2020

print(ensamble_name)
print(ssp)
    
#calibration
ys = y_alfa-6 #2011
ye = y_omega #2020
    
wspinup = 'wspinup'
#wspinup = ''

if cal == 'w5e5' and y_omega > 2018:
    print('w5e5 calibration out of bonds, last 2 years will be repeated')

spin_years = np.round(np.linspace(y_alfa-6, y_alfa, (y_alfa - (y_alfa-6)) * 12 + 1), 2) + 0.68 #0.67 stands for october
years = np.round(np.linspace(y_alfa, y_omega, (y_omega - y_alfa) * 12 + 1), 2) + 0.68 #0.67 stands for october
years = years[:-1] #get rid of the january y+1


# Out path
out_path = "/home/francesc/results/aneto/"
#resolution: n * n cells
#n = 5 #40 means 35m

# future
ys = y_ini = 2020
ys = 2020-6
ye = y_fin = 2099
y_ini=2020

#fcst years
years_fcst = np.round(np.linspace(y_ini, y_fin, (y_fin - y_ini) * 12 + 1), 2) + 0.68 #0.67 stands for october
