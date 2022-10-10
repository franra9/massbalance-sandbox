#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:18:00 2022

@author: francesc
"""
import numpy as np

# calibration period: 2011-2019 (w5e5)
y_alfa = 2011
y_omega = 2020
#cal = 'w5e5'
cal = 'isimip3b'

#'ssp126' #'ssp126' #"ssp245" does not exist!
ssp = 'ssp126' 
ssp = 'ssp370' 
ssp = 'ssp585'

if cal == 'w5e5': # do some assertion here
    ssp = ''
    
if cal == 'w5e5' and y_omega > 2019:
    Exception('w5e5 calibration out of bonds')

years = np.round(np.linspace(y_alfa, y_omega, (y_omega - y_alfa) * 12 + 1), 2) + 0.68 #0.67 stands for october
rho = 0.85 # 850kg/m3

# Out path
out_path = "/home/francesc/results/aneto/"
#resolution: n * n cells
n = 5