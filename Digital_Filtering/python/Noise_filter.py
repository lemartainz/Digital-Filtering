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


def filter_script(d_file, chan_num,
                  low_pass_a = 4e6, low_pass_b = 4e5,
                  high_pass_a = 2e6, high_pass_b = 2e5, 
                  calc_cut_numba_a = 5000., calc_cut_numba_b = 99999.99999999999):
    

    d = NA.digi_data(d_file, chan_num, convert_int = True)

    d.fft()
    
    d.set_low_pass(low_pass_a, low_pass_b)

    d.calc_low_pass()

    d.apply_lp_filter()

    d.set_high_pass(high_pass_a, high_pass_b)

    d.calc_high_pass()

    d.apply_hp_filter()

    d.calc_cut_numba(calc_cut_numba_a, calc_cut_numba_b)

    d.smooth_cuts(100.)

    d.apply_cuts()

    d.invert_corr()
    
    sl = d.get_t_slice(.1, .05)
    
    B.pl.figure()
    
    B.pl.plot(d.tall[sl], d.V[sl])

    B.pl.plot(d.tall[sl], np.real(d.V_c[sl]))

    B.pl.xlabel(r'Time ($\mu$s)')
    
    B.pl.ylabel('Voltage (V)')
    
    f_name = f'{d.name}_chan_num_{d.chan_num}_filtered.npz'
    np.savez(f_name, time = d.tall, signal = d.V_c)
    
#%% Directory Scraper

def scrape_data(dir_path):
    d_file = []
    for f in os.listdir(dir_path):
        if (fn.fnmatch(f, '*.hws') and fn.fnmatch(f, 'DAQ*')):
            d_file.append(f)
    print(d_file)
    return d_file
    
d_file = scrape_data('/Users/leo/Documents/GitHub/Digital-Filtering/Digital_Filtering/python')
chan_num = [0, 1, 2, 3]

#%% Filter Datasets
for i, j in enumerate(d_file):
    for k in chan_num:
        filter_script(j, k)
        #print(j,k)
