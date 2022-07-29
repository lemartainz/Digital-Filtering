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
                  low_pass_f = 3e6, low_pass_a = 3e5,
                  high_pass_f = 1e5, high_pass_a = 1e4, 
                  calc_cut_numba_a = 6000., calc_cut_numba_b = 99999.99999999999):
    

    d = NA.digi_data(d_file, chan_num, convert_int = True)

    d.fft()
    
    d.set_low_pass(low_pass_f, low_pass_a)

    d.calc_low_pass()

    d.apply_lp_filter()

    d.set_high_pass(high_pass_f, high_pass_a)

    d.calc_high_pass()

    d.apply_hp_filter()

    d.calc_cut_numba(calc_cut_numba_a, calc_cut_numba_b)

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

    
#%% Directory Scraper

def scrape_data(dir_path):
    d_file = []
    for f in os.listdir(dir_path):
        if (fn.fnmatch(f, '*.hws') and fn.fnmatch(f, 'DAQ*')):
            d_file.append(f)
    
    return d_file
    
d_file = scrape_data('G:\Github\Digital_Filtering\python')
chan_num = [0, 1, 2, 3]

#%% Filter Datasets
'''
for i, j in enumerate(d_file):
    for k in chan_num:
        filter_script(j, k)
        #print(j,k)
'''

#%% Optimize Parameters

low_pass_f = np.arange(1.1e6, 3.1e6, 1e5)
low_pass_a = np.arange(1e5, 3e5, 1e4)
high_pass_f = np.arange(0.5e5, 2.5e5, 1e4)
high_pass_a = np.arange(0.5e-2, 2.5e-2, 1e-3)


#%%
plot_raw(d_file[1], 0)
filter_script(d_file[1], 0, 3e6, 3e5, 1e5, 1e4, 6000., 99999.99999999999)
B.pl.xlim(0.116424, 0.116427)
B.pl.ylim(-0.15, 0.7)
B.pl.figure()


#%% HIGH PASS FREQ

for i in high_pass_f[:5]:
    filter_script(d_file[1], 0, 3e6, 3e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_f[0], high_pass_f[1], high_pass_f[2], high_pass_f[3], high_pass_f[4]], loc ='best')  
plot_raw(d_file[1], 0)
B.pl.figure()
    

for i in high_pass_a[5:10]:
    filter_script(d_file[1], 0, 3e6, 3e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_f[5], high_pass_f[6], high_pass_f[7], high_pass_f[8], high_pass_f[9]], loc ='best')
plot_raw(d_file[1], 0)
B.pl.figure()
    

for i in high_pass_f[10:15]:
    filter_script(d_file[1], 0, 3e6, 3e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_f[10], high_pass_f[11], high_pass_f[12], high_pass_f[13], high_pass_f[14]], loc = 'best')
plot_raw(d_file[1], 0)
B.pl.figure()


for i in high_pass_f[15:20]:
    filter_script(d_file[1], 0, 3e6, 3e5, i, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_f[15], high_pass_f[16], high_pass_f[17], high_pass_f[18], high_pass_f[19]], loc = 'best')
plot_raw(d_file[1], 0)   
B.pl.figure()


#%% HIGH PASS ALPHA


for i in high_pass_a[:5]:
    filter_script(d_file[1], 0, 5e6, 4e5, 1e5, i, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_a[0], high_pass_a[1], high_pass_a[2], high_pass_a[3], high_pass_a[4]], loc ='best')  
plot_raw(d_file[1], 0)
B.pl.figure()


for i in high_pass_a[5:10]:
    filter_script(d_file[1], 0, 5e6, 4e5, 1e5, i, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_a[5], high_pass_a[6], high_pass_a[7], high_pass_a[8], high_pass_a[9]], loc ='best')
plot_raw(d_file[1], 0)
B.pl.figure()
    
   
for i in high_pass_a[10:15]:
    filter_script(d_file[1], 0, 5e6, 4e5, 1e5, i, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_a[10], high_pass_a[11], high_pass_a[12], high_pass_a[13], high_pass_a[14]], loc = 'best')
plot_raw(d_file[1], 0)
B.pl.figure()


for i in high_pass_a[15:20]:
    filter_script(d_file[1], 0, 5e6, 4e5, 1e5, i, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([high_pass_a[15], high_pass_a[16], high_pass_a[17], high_pass_a[18], high_pass_a[19]], loc = 'best')
plot_raw(d_file[1], 0)                
B.pl.figure()

            
#%% LOW PASS ALPHA


for i in low_pass_a[:5]:
    filter_script(d_file[1], 0, 5e6, i, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_a[0], low_pass_a[1], low_pass_a[2], low_pass_a[3], low_pass_a[4]], loc ='best')  
plot_raw(d_file[1], 0)
B.pl.figure()


for i in low_pass_a[5:10]:
    filter_script(d_file[1], 0, 5e6, i, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_a[5], low_pass_a[6], low_pass_a[7], low_pass_a[8], low_pass_a[9]], loc ='best')
plot_raw(d_file[1], 0)
B.pl.figure()
    
    
for i in low_pass_a[10:15]:
    filter_script(d_file[1], 0, 5e6, i, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_a[10], low_pass_a[11], low_pass_a[12], low_pass_a[13], low_pass_a[14]], loc = 'best')
plot_raw(d_file[1], 0)
B.pl.figure()


for i in low_pass_a[15:20]:
    filter_script(d_file[1], 0, 5e6, i, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_a[15], low_pass_a[16], low_pass_a[17], low_pass_a[18], low_pass_a[19]], loc = 'best')
plot_raw(d_file[1], 0)   
B.pl.figure()
#%% LOW PASS FREQ  


for i in low_pass_f[:5]:
    filter_script(d_file[1], 0, i, 4e5, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_f[0], low_pass_f[1], low_pass_f[2], low_pass_f[3], low_pass_f[4]], loc ='best')  
plot_raw(d_file[1], 0)
B.pl.figure()
    

for i in low_pass_f[5:10]:
    filter_script(d_file[1], 0, i, 4e5, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_f[5], low_pass_f[6], low_pass_f[7], low_pass_f[8], low_pass_f[9]], loc ='best')
plot_raw(d_file[1], 0)
B.pl.figure()
    
   
for i in low_pass_f[10:15]:
    filter_script(d_file[1], 0, i, 4e5, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_f[10], low_pass_f[11], low_pass_f[12], low_pass_f[13], low_pass_f[14]], loc = 'best')
plot_raw(d_file[1], 0)
B.pl.figure()


for i in low_pass_f[15:20]:
    filter_script(d_file[1], 0, i, 4e5, 1e5, 1e4, 5000., 99999.99999999999)
    print(i)
    B.pl.xlim(0.116424, 0.116427)
    B.pl.ylim(-0.15, 0.6)
    B.pl.legend([low_pass_f[15], low_pass_f[16], low_pass_f[17], low_pass_f[18], low_pass_f[19]], loc = 'best')
plot_raw(d_file[1], 0)    
B.pl.figure()