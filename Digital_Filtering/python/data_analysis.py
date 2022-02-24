#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 08:06:16 2021

Testing digital filering of digitzer data

Start with simulated data

@author: boeglinw
"""

import numpy as np
import LT.box as B
import os
import h5py
import copy as C


from scipy.fft import fft, ifft, fftfreq, fftshift, ifftshift
from scipy.interpolate import UnivariateSpline

us= 1e-6

def set_file_info(filename):
     ddir, fname = os.path.split(filename)
     name, ext = os.path.splitext(fname)
     ddir += './'
     name = name
     ext = ext
     return ddir, name, ext


class digi_data:

    def __init__(self, file, chan_num, convert_int = False):
            self.dir, self.name, self.ext = set_file_info(file)
            print(("Open file : ",self.dir + self.name + self.ext))
            self.f = h5py.File(self.dir + self.name + self.ext, 'r')
            # get the data
            # ee information
            print("Get data")
            data_root = f'wfm_group0/traces/trace{chan_num:d}/'
            try:
                self.t0 = self.f[data_root + 'x-axis'].attrs['start']
                self.dt = self.f[data_root + 'x-axis'].attrs['increment']
                # measured data
                # scale dataset
                self.scale = self.f[data_root + 'y-axis/scale_coef'][()]
                # get the y dataset
                self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
            except:
                print("Problems loading data " + data_root)
                return
            self.ndata = self.nall
            self.ti = self.t0
            self.tf = self.t0 + (self.ndata-1)*self.dt
            self.tmin = self.t0
            self.tmax = self.t0 + (self.ndata-1)*self.dt
            if convert_int:
                 self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()].astype('int16')
            else:
                 self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()]

            # make the time axis
            print("Calculate data")
            self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
            self.V = self.scale[0] + self.scale[1]*self.ydata
            # select the data to be worked on
            self.f.close()
            print("Done")


#%% simulated data file

d_file = 'DAQ_160813-113140.hws'

d = digi_data(d_file, 0, convert_int = True)

t = d.tall
V = d.V

dt = np.diff(t[:2])[0]


N = 2* int( len(t)/2.)  # take only an even number of sample points

V_f = fft(V[:N])
V_i = ifft(V_f)

sp = fftshift(V_f)

f = fftshift(fftfreq(t[:N].shape[-1], d = dt))*us


#%% Smoothing with rolling average

# window size in nu ber of data points
kernel_size = 100
kernel = np.ones(kernel_size)/kernel_size


n_skip = 200000

# skip f = 0 part
sel = f!=0.


ps_s = np.convolve(np.abs(sp[sel][::n_skip]), kernel, mode = 'same')

f_s = f[sel][::n_skip]
ps_s_s = ps_s[::]

#
B.pl.plot(f_s, ps_s_s, '.')

# spline interpolation
sp_spl = UnivariateSpline(f_s, ps_s_s, k = 3)

ff = f[sel][::n_skip//5]
B.pl.plot(ff, sp_spl(ff))


#%% find large deviations
dsp = (sp - sp_spl(f))/sp_spl(f)
s_cut = (np.abs(f)>2.5)&(np.abs(dsp) > .5)

sp_c = C.copy(sp)
sp_c[sel & s_cut] = 0.

#%% fit a functional form

A = B.Parameter(500., 'A')
gam = B.Parameter(1.6, 'gam')

Ag = B.Parameter(5000., 'Ag')
sig = B.Parameter(1.6, 'sig')

Ag1 = B.Parameter(5000., 'Ag1')
sig1 = B.Parameter(1.6, 'sig1')

def gaus(x):
    y = Ag()*np.exp(-x**2/(2.*sig()**2))
    return y

def gaus1(x):
    y = Ag1()*np.exp(-x**2/(2.*sig1()**2))
    return y

def lorentz(x):
    g2 = 0.5*gam()
    y = A()/np.pi*g2/(x**2 + g2**2)
    return y

def signal(x):
    return gaus(x) + gaus1(x) + lorentz(x)

fit = B.genfit(signal,[Ag, sig, Ag1, sig1, A, gam], f[:-1][::n_skip] ,  ps_s)

#

B.pl.plot(f[:-1][::n_skip], fit.func(f[:-1][::n_skip]))

#%%
s_f = np.abs(sp[sel])/sp_spl(f[sel])

sp_c = np.zeros_like(sp)
sp_c[sel] = s_f*sp[sel]
sp_c[0] = sp[0]

V_c = ifft(ifftshift(sp_c))

#%%
istart = int(10e6); iend = int(5e5)
sl = slice(istart, istart + iend)

B.pl.plot(t[sl], np.real(V[sl]))


B.pl.plot(t[sl], np.real(V_i[sl]))

B.pl.plot(t[sl], np.real(V_c[sl]))