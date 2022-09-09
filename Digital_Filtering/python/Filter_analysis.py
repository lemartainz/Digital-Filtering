#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:41:34 2022

@author: leo
"""

import LT.box as B
import numpy as np
import Noise_algo as NA
import os
import fnmatch as fn


def filter_script(d_file, chan_num,
                  low_pass_f = 3e6, low_pass_a = 3e5,
                  high_pass_f = 1e5, high_pass_a = 1e4, 
                  calc_cut_numba_d = 6000., calc_cut_numba_p = 99999.99999999999):
    

    d = NA.digi_data(d_file, chan_num, convert_int = True)

    d.fft()
    
    d.set_low_pass(low_pass_f, low_pass_a)

    d.calc_low_pass()

    d.apply_lp_filter()

    d.set_high_pass(high_pass_f, high_pass_a)

    d.calc_high_pass()

    d.apply_hp_filter()

    d.calc_cut_numba(calc_cut_numba_d, calc_cut_numba_p)

    d.smooth_cuts(1000.)

    d.apply_cuts()

    d.invert_corr()
    
    sl = d.get_t_slice(.1, .05)

    B.pl.plot(d.tall[sl], np.real(d.V_c[sl]))
    
    B.pl.xlabel(r'Time (s)')
    
    B.pl.ylabel('Voltage (V)')
    
    f_name = f'{d.name}_chan_num_{d.chan_num}_filtered.npz'
    np.savez(f_name, time = d.tall, signal = d.V_c)
    print('Done with Filtering!')
    
    
def plot_raw(d_file, chan_num):
    
    d = NA.digi_data(d_file, chan_num, convert_int = True)
    
    sl = d.get_t_slice(.1, .05)
    
    B.pl.plot(d.tall[sl], d.V[sl], color = 'black')

    B.pl.xlabel(r'Time (s)')

    B.pl.ylabel('Voltage (V)')
    print('Done plotting raw data')

def scrape_data(dir_path):
    d_file = []
    for f in os.listdir(dir_path):
        if (fn.fnmatch(f, '*.hws') and fn.fnmatch(f, 'DAQ*')):
            d_file.append(f)
    
    return d_file
# windows directory    
#d_file = scrape_data('G:\Github\Digital_Filtering\python')
# mac directory
d_file = scrape_data('/Users/leo/Documents/GitHub/Digital-Filtering/Digital_Filtering/python')
chan_num = [0, 1, 2, 3]

#%% Filter Datasets


for i, j in enumerate(d_file):
    for k in chan_num:
        filter_script(j, k)
        
#%% Plot filtered

filt_data = B.get_file(input('Filtered data to plot: '))

#%%
u = input('INput: ')