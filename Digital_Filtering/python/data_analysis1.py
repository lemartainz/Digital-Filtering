#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 08:06:16 2021

Testing digital filering of digitzer data

Use real data

use rfft: fft of a real function

@author: boeglinw
"""

import numpy as np
import LT.box as B
import os
import h5py
import copy as C
from matplotlib.widgets import Cursor
import ffind_peaks as FP

# colored trajectories
import matplotlib.colors as COL
# colors list
color_table = COL.CSS4_COLORS
color_names = list(color_table.keys())

# color cycle
from itertools import cycle, islice
# select colors
color_selection = []


color_cycle = cycle(color_names)


from scipy.fft import rfft, irfft, rfftfreq
from scipy.interpolate import UnivariateSpline

us = 1e-6

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
            self.cuts = []
            self.first = True
            self.cid = None
            self.cid_k = None
            self.V_fc = None
            self.cut_fact = None
            self.lp_fact = None
            self.alpha = 1.
            self.fmax = 1e8

    def fft(self):
        t = self.tall
        V = self.V
        self.N = 2*(self.tall.size//2)
        self.V_f = rfft(V[:self.N])
        self.sp = np.abs(self.V_f)
        self.f = rfftfreq(t[:self.N].shape[-1], d = self.dt)
        self.df = np.diff(self.f[:2])[0]

    def plot_sp(self, skip):
        self.fig,self.ax  = B.pl.subplots()
        self.cursor = Cursor(self.ax, useblit=True, color='red', linewidth=1) # set special cursor
        self.ax.plot(self.f[::skip], self.sp[::skip])

    def plot_spc(self, skip):
        if self.V_fc is None:
            print('No corrected spectrum')
            return
        self.fig,self.ax  = B.pl.subplots()
        self.cursor = Cursor(self.ax, useblit=True, color='red', linewidth=1) # set special cursor
        self.ax.plot(self.f[::skip], np.abs(self.V_fc[::skip]))


    def reset_cuts(self):
        self.cuts = []


    def f_cut(self, skip):
        self.ml = None
        self.mh = None
        self.first = True
        self.xs = 0.
        self.xl = 0.
        self.plot_sp(skip)
        #self.cid_k = self.fig.canvas.mpl_connect('key_press_event', self.onkeypressed)
        #self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)


    def show_cuts(self):
        if self.fig is None:
            print('No figure, no plot !')
            return
        for x1, x2 in self.cuts:
            cc = next(color_cycle)
            ys,yl = B.pl.ylim()
            self.ax.plot([x1, x1], [0., 0.5*yl], color = cc)
            self.ax.plot([x2, x2], [0., 0.5*yl], color = cc)

    def add_cuts(self, method = 'M'):
        if self.fig is None:
            print('No figure !')
            return
        if method == 'K':
            self.method = method
            self.cid_k = self.fig.canvas.mpl_connect('key_press_event', self.onkeypressed)
        else:
            self.method = 'M'
            self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def done_cuts(self):
        # finished adding cuts
        if self.method == 'K':
            if self.cid_k is None:
                return
            self.fig.canvas.mpl_disconnect(self.cid_k)
        else:
            if self.cid is None:
                return
            self.fig.canvas.mpl_disconnect(self.cid)


    def onclick(self, event):
        self.xlim = B.pl.xlim()
        self.ylim = B.pl.ylim()        # entering cuts with mouse
        # double click add cut
        # right mouse button: finish
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
        if (event.button == 1):
            xp = event.xdata
            yp = event.ydata
            if self.first:
                if self.ml is not None:
                    self.ax.lines.remove(self.ml[0])
                self.ml = self.ax.plot([xp, xp], [0., yp], color = 'm')
                self.first = False
                self.xs = xp
            else:
                if self.mh is not None:
                    self.ax.lines.remove(self.mh[0])
                self.mh = self.ax.plot([xp, xp], [0., yp], color = 'r')
                self.xl = xp
                self.first = True
            B.pl.xlim(self.xlim)
            B.pl.ylim(self.ylim)
            B.pl.draw()
            print(f'xs = {self.xs}, xl = {self.xl}')
        elif (event.button == 3):
            # add cut
            if self.xs > self.xl:
                # make sure xs < xl
                xx = self.xs
                self.xs = self.xl
                self.xl = xx
            print(f'add cut = ({self.xs},{self.xl})')
            self.cuts.append((self.xs, self.xl))



    def onkeypressed(self, event):
        # entering cuts with keys
        print(f'Key pressed = {event.key}, xdata = {event.xdata}, ydata = {event.ydata}')
        self.xlim = B.pl.xlim()
        self.ylim = B.pl.ylim()
        yh = self.ylim[-1]
        if event.key == ' ':
            # space bar left marker
            if self.ml is not None:
                self.ax.lines.remove(self.ml[0])
            xs = event.xdata
            self.ml = self.ax.plot([xs, xs], [0., 0.5*yh], color = 'r')
            B.pl.xlim(self.xlim)
            B.pl.ylim(self.ylim)
            B.pl.draw()
            print(f'xs = {xs}')
            self.xs = xs
        elif event.key == 'ctrl+ ':
            if self.mh is not None:
                self.ax.lines.remove(self.mh[0])
            xl = event.xdata
            self.mh = self.ax.plot([xl, xl], [0., 0.5*yh], color = 'r')
            B.pl.xlim(self.xlim)
            B.pl.ylim(self.ylim)
            B.pl.draw()
            print(f'xl = {xl}')
            self.xl= xl
        elif event.key == 'enter':
            # add cut
            if self.xs > self.xl:
                # make sure xs < xl
                xx = self.xs
                self.xs = self.xl
                self.xl = xx
            print(f'add cut = ({self.xs},{self.xl})')
            self.cuts.append((self.xs, self.xl))
        elif event.key == 'ctrl+enter':
            print('Done entering cuts')
            self.fig.canvas.mpl_disconnect(self.cid_k)
        else:
            return

    def calc_cuts(self):
        # calculate the cut factor to be applied later
        self.cut_fact = np.ones_like(self.f)
        for xs, xl in self.cuts:
            self.cut_fact [(xs <self.f)&(self.f<xl)] = 0.

    def smooth_cuts(self, sigma):
        # smooth the cut edges
        kernel_size = int(5.*sigma/self.df)
        x = np.linspace(-5., 5., kernel_size)
        kernel0 = np.exp(-x**2/2.)
        kernel = kernel0/kernel0.sum()
        self.cut_fact =  np.convolve(self.cut_fact, kernel, mode = 'same')


    def set_low_pass(self, fmax, alpha):
        self.fmax = fmax
        self.alpha = alpha


    def calc_low_pass(self):
        # pass freq. up to self.xmax the steepness of the cut is controlled by alpha
        # alpha width in Hz for step
        # index array of values to cut
        self.lp_fact = 0.5*(1. - np.tanh( 1./self.alpha * (self.f - self.fmax) ))

    def apply_cuts(self):
        if self.cut_fact is None:
            print('cuts not calculated (use calc_cuts)')
            return
        if self.V_fc is None:
            self.V_fc = C.copy(self.V_f)
        self.V_fc *= self.cut_fact

    def reset_corr(self):
        self.V_fc = C.copy(self.V_f)

    def invert_corr(self):
        if self.V_fc is None:
            return
        self.V_c = irfft(self.V_fc)

    def apply_lp_filer(self):
        if self.lp_fact is None:
            print('low pass filer not calculated (use low_pass)')
            return
        if self.V_fc is None:
            self.V_fc = C.copy(self.V_f)
        self.V_fc *= self.lp_fact

    def save_filters(self, f_filename):
        o = open(f_filename, 'w')
        o.write(f'#\ alpha = {self.alpha} \n')
        o.write(f'#\ fmax = {self.fmax} \n')
        o.write('#! xs[f,0]/ xl[f,1]/ \n')
        for xs, xl in self.cuts:
            o.write(f'{xs} {xl} \n')
        o.close()
 
    def load_filters(self, f_filename):
         fd = B.get_file(f_filename)
         self.alpha = fd.par['alpha']
         self.fmax = fd.par['fmax']
         self.cuts = list(zip(fd['xs'],fd['xl']))

    def get_t_slice(self, t_start, delta_t):
        istart = int(t_start/self.dt); iend = int(delta_t/self.dt)
        return slice(istart, istart + iend)

    def get_f_slice(self, f_start, delta_f):
        istart = int(f_start/self.df); iend = int(delta_f/self.df)
        return slice(istart, istart + iend)
    
    def find_peaks(yval, ystep, xval = None, \
               xmin = None, xmax = None, limits = None, \
               power = 5,\
               get_window = None ):
     # find peaks in an array of data
     # get_window is a function that returns the slice of data between xmin and xmax
     # peak finding
     nmin = 0
     nmax = 0
     MAXTAB=[]
     MINTAB=[]
     if ((xmin == None) and (xmax == None)) or (not limits) or (get_window == None):
          print("No limits present, analyze all data")
          results = np.zeros((2,), dtype = 'int32')
          pmin = np.zeros((len(yval)//5, ), dtype='int32')
          pmax = np.zeros((len(yval)//5, ), dtype='int32')
          try:
              FP.find_peaks(len(yval), ystep, yval, results, pmin, pmax)
          except:
              print("problem with peak finding")
              return []
          nmin = results[0]
          nmax = results[1]
     else:
          # get the window
          print(("Analyze data between ", xmin, " and ", xmax))
          sl = get_window( xmin, xval, xmax)
          # progress dialog
          results = np.zeros((2,), dtype = 'int32')
          pmin = np.zeros((len(yval[sl])//5, ), dtype='int32')
          pmax = np.zeros((len(yval[sl])//5, ), dtype='int32')
          try:
              FP.find_peaks(len(yval[sl]), ystep, yval[sl], results, pmin, pmax)
          except:
               print("problem with peak finding")
               return []
          nmin = results[0]
          nmax = results[1]
     MAXTAB.append( xval[pmax[:nmax]] )
     MAXTAB.append( yval[pmax[:nmax]] )
     MAXTAB.append(pmax[:nmax])
     MINTAB.append( xval[pmin[:nmin]] )
     MINTAB.append( yval[pmin[:nmin]] )
     MINTAB.append(pmin[:nmin])
     return [MAXTAB,MINTAB]


#%% simulated data file

d_file = 'DAQ_160813-113140.hws'
#d_file = '../MAST_data/Aug_13/DAQ_190813-112521.hws'

d = digi_data(d_file, 0, convert_int = True)

d.fft()

#%% load filter and apply
# ready to apply cuts etc.

d.load_filters('test_filer.data')
d.calc_cuts()
d.smooth_cuts(250.)
d.calc_low_pass()
d.apply_lp_filer()
d.apply_cuts()
d.invert_corr()


#%%
sl = d.get_t_slice(.1, .05)

#%%
B.pl.plot(d.tall[sl], d.V[sl])

B.pl.plot(d.tall[sl], np.real(d.V_c[sl]))

#%% save corrected data

f_name = f'{d.name}_filtered.npz'
np.savez(f_name, time = d.tall, signal = d.V_c)
