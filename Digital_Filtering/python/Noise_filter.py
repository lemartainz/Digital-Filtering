#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:56:13 2022

@author: leo
"""

import LT.box as B
import numpy as np
import Noise_algo as NA
import os
import fnmatch as fn
import tqdm 

def filter_script(d_file, chan_num,
                  low_pass_a = 4e6, low_pass_f = 4e5,
                  high_pass_a = 1e5, high_pass_f = 1e4, 
                  calc_cut_numba_a = 5000., calc_cut_numba_b = 99999.99999999999):
    

    d = NA.digi_data(d_file, chan_num, convert_int = True)

    d.fft()
    
    d.set_low_pass(low_pass_a, low_pass_f)

    d.calc_low_pass()

    d.apply_lp_filter()

    d.set_high_pass(high_pass_a, high_pass_f)

    d.calc_high_pass()

    d.apply_hp_filter()

    d.calc_cut_numba(calc_cut_numba_a, calc_cut_numba_b)

    d.smooth_cuts(100.)

    d.apply_cuts()

    d.invert_corr()
    
    sl = d.get_t_slice(.1, .05)
    
    # B.pl.figure()
    
    B.pl.plot(d.tall[sl], d.V[sl])

    B.pl.plot(d.tall[sl], np.real(d.V_c[sl]))

    B.pl.xlabel(r'Time ($\mu$s)')
    
    B.pl.ylabel('Voltage (V)')
    
    f_name = f'{d.name}_chan_num_{d.chan_num}_filtered.npz'
    np.savez(f_name, time = d.tall, signal = d.V_c)
    print('Done with Filtering!')
    
#%% Directory Scraper

def scrape_data(dir_path):
    d_file = []
    for f in os.listdir(dir_path):
        if (fn.fnmatch(f, '*.hws') and fn.fnmatch(f, 'DAQ*')):
            d_file.append(f)
    
    return d_file
    
d_file = scrape_data('/Users/leo/Documents/GitHub/Digital-Filtering/Digital_Filtering/python')
chan_num = [0, 1, 2, 3]

#%% Filter Datasets
'''
for i, j in enumerate(d_file):
    for k in chan_num:
        filter_script(j, k)
        #print(j,k)
'''

#%% Optimize Parameters

low_pass_a = np.arange(3e6, 5e6, 1e5)
low_pass_f = np.arange(3e5, 5e5, 1e4)
high_pass_a = np.arange(0.5e5, 2.5e5, 1e4)
high_pass_f = np.arange(0.5e4, 2.5e4, 1e3)

filter_script(d_file[0], 0, 5e6, 4e5, 1e5, 1e4, 5000., 99999.99999999999)
B.pl.figure()
for i in high_pass_a[:5]:
    filter_script(d_file[0], 0, 5e6, 4e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    p = []
    p.append(i)
    B.pl.legend([50000, 60000, 70000, 80000, 90000])  
B.pl.figure()
    
for i in high_pass_a[5:10]:
    filter_script(d_file[0], 0, 5e6, 4e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    p = []
    p.append(i)
    B.pl.legend([100000, 110000, 120000, 130000, 140000])
B.pl.figure()
    
    
for i in high_pass_a[10:15]:
    filter_script(d_file[0], 0, 5e6, 4e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    p = []
    p.append(i)
    B.pl.legend([150000, 160000, 170000, 180000, 190000])
B.pl.figure()

for i in high_pass_a[15:20]:
    filter_script(d_file[0], 0, 5e6, 4e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    p = []
    p.append(i)
    B.pl.legend([200000, 210000, 220000, 230000, 240000])