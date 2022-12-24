#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:18:00 2022
# calibration params file
@author: francesc
"""
import sys

#print('Number of arguments:'), print(len(sys.argv), 'arguments.')
#print('Argument List:'), print(str(sys.argv))

# thing to run many simulations at a time
#ens = int(sys.argv[1])
#ss = int(sys.argv[2])
#n = int(sys.argv[3])

import numpy as np

# calibration period: 2011-2019 (w5e5)

cal = 'w5e5'
#cal = 'isimip3b'

ensamble_names = ['ukesm1-0-ll_r1i1p1f2', 'gfdl-esm4_r1i1p1f1', 'ipsl-cm6a-lr_r1i1p1f1',
       'mpi-esm1-2-hr_r1i1p1f1', 'mri-esm2-0_r1i1p1f1', 'w5e5']

ens=5 #2
ss=1
n=40

ensamble_name=ensamble_names[ens]
#ensamble_name='gfdl-esm4_r1i1p1f1'

#'ssp126' #'ssp126' #"ssp245" does not exist!
ssps = ['ssp126', 'ssp370', 'ssp585']
ssp = ssps[ss]
#ssp = 'ssp370' 
#ssp = 'ssp585'
# paralelization:
ncores = 4



if cal == 'w5e5': # do some assertion here
    ssp = ''
    y_alfa = 2011
    y_omega = 2018
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
rho = 0.85 # 850kg/m3

# Out path
out_path = "/home/francesc/results/aneto/"
