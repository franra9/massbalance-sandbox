#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 19:10:15 2022

@author: francesc
"""
import numpy as np
import pandas as pd
import params as params
import os

data_path='/home/francesc/data/aneto_glacier/climate/AWS_aneto/'
d11=pd.read_csv(f'{data_path}/Datos20110.csv', encoding = "ISO-8859-1")
d12=pd.read_csv(f'{data_path}/Datos20120.csv', encoding = "ISO-8859-1")
d13=pd.read_csv(f'{data_path}/Datos20130.csv', encoding = "ISO-8859-1")
d14=pd.read_csv(f'{data_path}/Datos20140.csv', encoding = "ISO-8859-1")
d15=pd.read_csv(f'{data_path}/Datos20150.csv', encoding = "ISO-8859-1")
d16=pd.read_csv(f'{data_path}/Datos20160.csv', encoding = "ISO-8859-1")
d17=pd.read_csv(f'{data_path}/Datos20170.csv', encoding = "ISO-8859-1")
d18=pd.read_csv(f'{data_path}/Datos20180.csv', encoding = "ISO-8859-1")
d19=pd.read_csv(f'{data_path}/Datos20190.csv', encoding = "ISO-8859-1")
#d20=pd.read_csv(f'{data_path}/Datos20200.csv', encoding = "ISO-8859-1")

pd.read_table()

m1=d11['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m2=d12['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m3=d13['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m4=d13['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
#d14['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m6=d15['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m7=d16['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m8=d17['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m9=d18['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()
m10=d19['Temperatura (§C)'].str.replace(r" ", "nan").astype(float).mean()

np.mean([m1,m2,m3,m4,m6,m7,m8,m9,m10])

#-1.2970044627448898

from climate import get_climate

years=np.linpace(2011,2020,72)

clim=get_climate(years,3050)

dum=[]
for i in np.arange(0, len(clim.index)):
    dum.append(np.mean(clim.temp.iloc[i]))

print(np.mean(dum))
#0.4128301154469802







