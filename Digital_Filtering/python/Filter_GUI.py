# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 12:28:26 2022

@author: mleon
"""
from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np
import LT.box as B
from LT.parameterfile import pfile
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as pl
from matplotlib.figure import Figure
import os
import sys
import os.path
import Noise_algo as NA
import h5py




#------------------------------------------------------------------------------
# Useful functions
#------------------------------------------------------------------------------

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



class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates=True)
        self.rescale_data = None
        self.N = 1000

        self.Ncontrol= QtWidgets.QLineEdit(self)
        self.Ncontrol.setValidator(QtGui.QIntValidator(self.Ncontrol))
        self.Ncontrol.setFixedWidth(50)
        self.Ncontrol.setText(str(self.N))

        self.refr=QtWidgets.QPushButton(QtGui.QIcon('refresh.png'), None, self)

        self.refr.clicked.connect(self.getN)

        #self.Ncontrol.textChanged.connect(self.getN)

        self.Nlabel=QtWidgets.QLabel('Data points on Fig, N=', self)


        self.Ncact=self.addWidget(self.Nlabel)
        self.Nlact=self.addWidget(self.Ncontrol)
        self.Nrefr=self.addWidget(self.refr)

    # to remove matplotlib toolbar status line
    def set_message(self, msg):
        pass

    def getN(self):
      self.N = int(self.Ncontrol.text())
      try:
          self.plot_rescaled()
      except:
          print('No plot yet to rescale')
    def back(self, *args):
        NavigationToolbar2QT.back(self, *args)
        self.plot_rescaled()

    def forward(self, *args):
        NavigationToolbar2QT.forward(self, *args)
        self.plot_rescaled()

    def release_pan(self, event):
        NavigationToolbar2QT.release_pan(self, event)
        self.plot_rescaled()

    def home(self, *args):
        NavigationToolbar2QT.home(self, *args)
        self.plot_rescaled()

    def release_zoom(self, event):
        NavigationToolbar2QT.release_zoom(self, event)
        self.plot_rescaled()

    #rescale and replot
    def plot_rescaled(self):
        # get the parent window
        parent = self.canvas.parent()
        #
        if self.rescale_data is None:
            # No data to plot
            return
        t=self.rescale_data[0]
        V=self.rescale_data[1]
        # get the axes of the parent window
        ax=parent.axes
        (xmin, xmax)=ax.get_xlim()
        rng=xmax-xmin
        #get original data which is in plotting region
        org_points=np.where(((xmin - 2*rng) < t) & (t < (xmax + 2*rng)))
        if (xmin- 2*rng)<t[0] or (xmax+2*rng)>t[-1]:
            k= (-max(xmin-2*rng, t[0])+min(t[-1], xmax+2*rng))/(xmax-xmin)
        else: k=5
        to=t[org_points]
        Vo=V[org_points]
        n=1 #take every nth point from data
        if len(to)/k> self.N:
            n=int(len(to)/k/self.N)
        tcut=to[::n]
        Vcut=Vo[::n]

        ax.set_prop_cycle(None)
        ax.lines.remove(parent.all_plot[0])
        ax.autoscale(enable=False, axis='x', tight=True)
        ax.set_autoscaley_on(False)

        if parent.par['draw_lines']:
            parent.all_plot=ax.plot(tcut,Vcut, zorder=0)
        else:
            parent.all_plot=ax.plot(tcut,Vcut, '.', zorder=0)
        ax.figure.canvas.draw()
        
        
class NumberDialog(QtWidgets.QDialog):
     def __init__(self, data, parent = None, title = 'Enter parameters', labels=['First Label','Second Label'],
                  keys=['first','second'], about_txt = 'about_txt'):

          self.parent=parent
          self.keys=keys
          self.data=data
          QtWidgets.QDialog.__init__(self, parent)
          self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
          self.setWindowTitle(title)
          self.layout=QtWidgets.QGridLayout(self)
          about = QtWidgets.QLabel(about_txt, self)
          ok=QtWidgets.QPushButton('Ok', self)
          ok.setDefault(True)
          cancel=QtWidgets.QPushButton('Cancel', self)
          ok.clicked.connect(self.OnOk)
          cancel.clicked.connect(self.OnCancel)
          sepline=QtWidgets.QFrame()
          sepline.setFrameShape(QtWidgets.QFrame.HLine)
          sepline.setFrameShadow(QtWidgets.QFrame.Sunken)

          self.layout.addWidget(about, 0, 0)
          self.layout.addWidget(sepline, 1, 0, 1, 2)


          # loop over keys to add controls and validators
          nrow = len(keys)+1
          #qle - qlinedit dictionary to retrieve data later
          self.qle={}
          for i, key in enumerate(keys):
               Val_l  = QtWidgets.QLabel(labels[i], self)
               Val_t  = QtWidgets.QLineEdit(self)
               Val_t.setValidator(QtGui.QDoubleValidator(self))

               Val_t.textChanged.connect(self.check_state)
               Val_t.setText(data.get(key))
               self.layout.addWidget(Val_l, i+2, 0)
               self.layout.addWidget(Val_t, i+2, 1)

               self.qle[key]=Val_t

          self.layout.addWidget(ok, nrow+2, 0)
          self.layout.addWidget(cancel, nrow+2, 1)

     def check_state(self, *args, **kwargs):

        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = '#c4df9b' # green
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)

     def OnOk(self):
         for i, key in enumerate(self.keys):
             try:
                 self.data[key]=self.qle[key].text()
             except:
                 mb=QtWidgets.QMessageBox(self)
                 mb.setWindowTitle('Entry error')
                 mb.setText("Please enter the missing parameter")
                 mb.exec_()
                 self.qle[key].setFocus()
                 return
         self.close()

     def OnCancel(self):
         self.destroy()



#------------------------------------------------------------------------------
# Main Plot Frame for raw and filtered data
#------------------------------------------------------------------------------



class Main_PlotFrame(QtWidgets.QMainWindow):
    
    def __init__(self, parent):
        QtWidgets.QMainWindow.__init__(self)
        
        self.setWindowTitle('Digital Filter')
        
        # parameters
        
        self.par = {}
        self.par['chan_num'] = 0
        self.par['low_pass_f'] = 3e6
        self.par['low_pass_a'] = 3e5
        self.par['high_pass_f'] = 1e5
        self.par['jigh_pass_a'] = 1e4
        self.par['calc_cut_numba_d'] = 6000.
        self.par['calc_cut_numba_p'] = 99999.99999999999
        self.par['convert_int'] = True
        self.par["draw_lines"] = False
        
        # Change background color
        
        self.palette=QtGui.QPalette()
        self.palette.setColor(QtGui.QPalette.Background,QtCore.Qt.white)
        
        # menu bar
        self.mbar = self.menuBar()
        
        # setup menues

        # create status bar
        self.stBar1 = QtWidgets.QLabel('No file loaded')
        self.stBar1.setFrameStyle(2)
        self.stBar1.setFrameShadow(48)
        self.stBar2 = QtWidgets.QLabel('No plot yet')
        self.stBar2.setFrameStyle(2)
        self.stBar2.setFrameShadow(48)

        self.statusBar().addWidget(self.stBar1, 1)
        self.statusBar().addWidget(self.stBar2, 1)

        # setup figure
        self.figure = Figure(figsize=(10,6))
        self.axes = self.figure.add_subplot(111)
        self.figure_canvas = FigureCanvas(self.figure)
        
        # inital settings
        self.ndata = 0
        self.t = None
        self.V = None
        self.Vinv = None
        self.t_slice = None
        self.V_slice = None
        self.hwsproc = [] #hws files selection
        self.hwsnew = []
        self.hwstodo = []
        self.datadir = './'
        #load list of processed files
        self.proclistload()

        # file information data
        self.parent = parent
        self.full_name = None
        self.dir = None
        self.par_dir = None
        self.hist_dir = None
        self.res_dir = None
        self.name = None
        self.ext = None
        self.fileDataOK = False  
    
    
        self.setCentralWidget(self.figure_canvas)    
        # toolbar
        self.toolbar = NavigationToolbar(self.figure_canvas, self)
        self.addToolBar(self.toolbar)


        self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      drawtype='box', useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
    
        self.setPalette(self.palette)
        self.setGeometry(150,150,800,600)
        self.createMenubar()
        self.show()
    def my_plot(self, *args, **kwargs):

        N=kwargs.pop('N', self.toolbar.N)
        ax=self.axes

        #saves input data into toolbar class object and reuses later
        self.toolbar.rescale_data = args

        #dropping points for first plot
        t=args[0]
        V=args[1]
        rest=args[2:]
        n=1
        #take every nth point from data
        if len(t)> N:
            n=int(len(t)/N)
        tcut=t[::n]
        Vcut=V[::n]

        ax.autoscale(enable=True, axis='x', tight=True)
        return ax.plot(tcut, Vcut, *rest, **kwargs)

    def proclistload(self): #load list of processed files
        path_to_watch = self.datadir
        dircont = os.listdir (path_to_watch)

        if "ProcFiles.data" in dircont:
            print("\nLoad list of already processed files.")
            fopen = open(self.datadir+"//ProcFiles.data")
            self.hwsproc=[x.strip('\n') for x in fopen.readlines()]

        else:
            self.hwsproc=[]
            print("\nNo files were processed yet in current data directory.")

    def proclistsave(self): #save list of processed files
        path_to_watch = self.datadir
        dircont = os.listdir (path_to_watch)

        if "ProcFiles.data" in dircont:
            print("\nAppend list of processed files.")
            fopen = open(self.datadir+"//ProcFiles.data", 'a')
            for item in self.hwsproc:
                fopen.write("%s\n" % item)
            fopen.close()
        else:
            fopen = open(self.datadir+"//ProcFiles.data", 'w')
            for item in self.hwsproc:
                fopen.write("%s\n" % item)
            fopen.close()
            print("\nCreate list of processed files.")


    def menuData(self): # data for the menu
          return(
               ("&File",                                        # menu label
                ("&Open", "Select data files", self.OnSelectFiles),  # menu items consisting of: label, text, handler
                ("&Reload", "Reload data files", self.OnReload),
                ("&Load Parameters", "Load parameters", self.OnLoadParameters),
                ("&Save Parameters", "Save parameters", self.OnSaveParameters),
                ("&Quit", "Quit program", self.OnCloseWindow)), # label, status and handler
               #
               ("&Actions",
                ("&Plot", "Plot data", self.OnPlot),
                (None, None, None),  # creates a separator bar in the menu
                ("&Clear Figure", "Clear figure", self.OnClear)),
               #
               ("&Parameters",
                ("&Data directory", "Set data directory", self.OnSetScanDir),
                (None, None, None),
                ("Detector Channel", "Set detector channel", self.OnSelectChannel),
                (None, None, None))
                )

    def createMenubar(self):
         group=QtWidgets.QActionGroup(self)
         for eachMenuData in self.menuData():

              menuLabel = eachMenuData[0]
              menuItems = eachMenuData[1:]


              menu = QtWidgets.QMenu(menuLabel, self)
              for eachItem in menuItems:
                   if not eachItem[0]:
                       menu.addSeparator()
                       continue
                   subMenu=QtWidgets.QAction(eachItem[0],self)
                   if len(eachItem) == 4: #never true
                        if eachItem[3]=='CHECK' or eachItem[3]=='RADIO':
                            subMenu.setCheckable(True)
                        if eachItem[3]=='RADIO':
                            subMenu.setActionGroup(group)


                   if 'Measure' in eachItem[0]: subMenu.setChecked(self.par['measure'])
                   if 'Convert Integer' in eachItem[0]: subMenu.setChecked(self.par['convert_int'])
                   if 'SetLimits' in eachItem[0]: subMenu.setChecked(self.par['limits'])
                   if 'Use Limits' in eachItem[0]: subMenu.setChecked(self.par['use_limits'])
                   if 'Plot Lines' in eachItem[0]: subMenu.setChecked(self.par['draw_lines'])
                   if 'Auto Histogram' in eachItem[0]: subMenu.setChecked(self.par['auto_histo'])
                   if 'Histogram Points' in eachItem[0]: subMenu.setChecked(self.par['plot_histo_points'])
                   if 'Filtered' in eachItem[0]: subMenu.setChecked(self.par['filtered'])


                   subMenu.triggered.connect(eachItem[2])
                   menu.addAction(subMenu)



              self.mbar.addMenu(menu)





    def set_file_info(self,filename):
          dir, fname = os.path.split(filename)
          name, ext = os.path.splitext(fname)
          self.dir = dir+'//'
          self.name = name
          self.ext = ext
          # that's it
          
    def OnSelectFiles(self):
         # Create a file-open dialog in the current directory
         filetypes = '*.hws'
         # use a regular file dialog
         if self.dir == None:
             self.dir = os.getcwd()

         file_dlg=QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file', self.dir, filetypes )

         if file_dlg[0] != '':
             # User has selected something, get the path, set the window's title to the path
             filename=file_dlg[0]
             # store relevant file information
             self.set_file_info(filename)
         else:
             print("so, you changed your mind, I will do nothing")
             filename = None
             return


         # get channel number
         chan_num = "%0d"%(int(self.par["chan_num"]))
         self.stBar1.setText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num))
         # open file
         # open files
         print(("Open file : ",self.dir + self.name + self.ext))
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # ee information
         print("Get data")
         data_root = 'wfm_group0/traces/trace' + chan_num + '/'
         try:
             self.t0 = self.f[data_root + 'x-axis'].attrs['start']
             self.dt = self.f[data_root + 'x-axis'].attrs['increment']
             # measured data
             # scale dataset
             self.scale = self.f[data_root + 'y-axis/scale_coef'][()]
             # get the y dataset
             self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
         except:
             mb=QtWidgets.QMessageBox(self)
             mb.setText("Problems loading data " + data_root)
             mb.exec_()
             return
         self.ndata = self.nall
         self.ti = self.t0
         self.tf = self.t0 + (self.ndata-1)*self.dt
         self.par["tmin"] = self.t0
         self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         # WB temp solution 8/13/13
         # set the type as 16 bit signed, this is not good the data type should come from
         # the hdf file
         if self.par["convert_int"]:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()].astype('int16')
         else:
             self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()]

         # make the time axis
         print("Calculate data")
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print("Select data")
         self.select_data()
         print("Done")
             
         # else:
         #     # load npz data file
         #     print(("Open npz data file : ",self.dir + self.name + self.ext))
         #     self.f = np.load(self.dir + self.name + self.ext)
         #     d = self.f
         #     print("Get npz data")
         #     self.t0 = d['time'][0]
         #     self.dt = np.diff(d['time'])[0]
         #     self.scale = [0.,1.]
         #     self.nall = len(d['time'])
         #     self.ndata = self.nall
         #     self.ti = self.t0
         #     self.tf = self.t0 + (self.ndata-1)*self.dt
         #     self.par["tmin"] = self.t0
         #     self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         #     self.tall = d['time']
         #     self.ydata = d['signal']
         #     print("Select data")
         #     self.select_data()
         #     print("Done")

    def OpenFile(self, fpname):

         self.set_file_info(fpname)
         chan_num = "%0d"%(int(self.par["Detector channel"]))
         self.stBar1.setText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num))
         # open files
         
         # real data
         print(("Open file : ",self.dir + self.name + self.ext))
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # time information
         print("Get data")
         data_root = 'wfm_group0/traces/trace' + chan_num + '/'
         try:
             self.t0 = self.f[data_root + 'x-axis'].attrs['start']
             self.dt = self.f[data_root + 'x-axis'].attrs['increment']
             # measured data
             # scale dataset
             self.scale = self.f[data_root + 'y-axis/scale_coef'][()]
             # get the y dataset
             self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
         except:
             mb=QtWidgets.QMessageBox(self)
             mb.setText("Problems loading data " + data_root)
             mb.exec_()
             return
         self.ndata = self.nall
         self.ti = self.t0
         self.tf = self.t0 + (self.ndata-1)*self.dt
         self.par["tmin"] = self.t0
         self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         
         if self.par["convert_int"]:
                self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()].astype('int16')
         else:
                self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()]

         # make the time axis
         print("Calculate data")
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print("Select data")
         self.select_data()
         print("Done")
         # else:
         #     # load npz data file
         #     print(("Open npz data file : ",self.dir + self.name + self.ext))
         #     self.f = np.load(self.dir + self.name + self.ext)
         #     d = self.f
         #     print("Get npz data")
         #     self.t0 = d['time'][0]
         #     self.dt = np.diff(d['time'])[0]
         #     self.scale =  [0.,1.]
         #     self.nall = len(d['time'])
         #     self.ndata = self.nall
         #     self.ti = self.t0
         #     self.tf = self.t0 + (self.ndata-1)*self.dt
         #     self.par["tmin"] = self.t0
         #     self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         #     self.tall = d['time']
         #     self.ydata = d['signal']
         #     print("Select data")
         #     self.select_data()
         #     print("Done")

    def OnReload(self):
         # reload data file
         # get channel number
         chan_num = "%0d"%(int(self.par["Detector channel"]))
         self.stBar1.setText('Current file : %s / Channel : %s'%(self.name+self.ext, chan_num))
         # open file
         print(("Open file : ",self.dir + self.name + self.ext))
         self.f = h5py.File(self.dir + self.name + self.ext, 'r')
         # get the data
         # time information
         print("Get data")
         data_root = 'wfm_group0/traces/trace' + chan_num + '/'
         try:
             self.t0 = self.f[data_root + 'x-axis'].attrs['start']
             self.dt = self.f[data_root + 'x-axis'].attrs['increment']
             # measured data
             # scale dataset
             self.scale = self.f[data_root + 'y-axis/scale_coef'][()]
             # get the y dataset
             self.nall = self.f[data_root + 'y-axis/data_vector/data'].shape[0]
         except:
             mb=QtWidgets.QMessageBox(self)
             mb.setText("Problems loading data " + data_root)
             mb.exec_()
             return
         self.ndata = self.nall
         self.ti = self.t0
         self.tf = self.t0 + (self.ndata-1)*self.dt
         # self.par["tmax"] = self.t0 + (self.ndata-1)*self.dt
         # WB temp solution 8/13/13
         # set the type as 16 bit signed, this is not good the data type should come from
         # the hdf file
         if self.par["convert_int"]:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()].astype('int16')
         else:
              self.ydata = self.f[data_root + 'y-axis/data_vector/data'][()]
         # make the time axis
         print("Calculate data")
         self.tall = self.t0 + self.dt*np.arange(self.ydata.shape[0], dtype = float)
         # select the data to be worked on
         self.f.close()
         print("Select data")
         self.select_data()
         print("Done")

    def OnLoadParameters(self):
          # Create a file-open dialog in the current directory
          filetypes = '*.data'
          # use a regular file dialog
          file_dlg=QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file', self.dir, filetypes )
          if file_dlg[0] != '':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg[0]
               # store relevant file information
               print(("will open " , filename))
          else:
               print("so, you changed your mind, I will do nothing")
               filename = None
               return
          p = pfile(filename)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print((key, " ", self.par[key]))
          # now set the check marks
          # get the menus

          menus = [i.menu() for i in self.mbar.actions()]
          items = None
          item_dict = {}
          # search for the options menu and get the contents
          for m in menus:
               if m.title() == "&Options":
                    items = m.actions()
                    for item in items:
                         if item.text() == "":
                              continue
                         else:
                              key = item.text().strip()[1:]
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["Convert Integer"].setChecked(self.par["convert_int"]); print(("Convert Integer = ", self.par["convert_int"]))
               item_dict["SetLimits"].setChecked(self.par["limits"]); print(("SetLimits = ", self.par["limits"]))
               item_dict["Measure"].setChecked(self.par["measure"]); print(("Measure = ", self.par["measure"]))
               item_dict["Use Limits"].setChecked(self.par["use_limits"]); print(("Use Limits = ",self.par["use_limits"]))
               item_dict["Auto Histogram Limits"].setChecked(self.par["auto_histo"]); print(("Auto Histogram Limits = ", self.par["auto_histo"]))
               item_dict['Histogram Points'].setChecked(self.par['plot_histo_points']);print(("Plot Histogram Points = ", self.par["plot_histo_points"]))
               item_dict["Filtered"].setChecked(self.par["filtered"]); print(("Filtered = ", self.par["filtered"]))
               item_dict["Plot Lines"].setChecked(self.par["draw_lines"]); print(("draw_lines = ", self.par["draw_lines"]))
          self.select_data()
          # all done

    def LoadParameters(self,fpname):

          p = pfile(fpname)
          for key in self.par:
               self.par[key] = p.get_value(key)
               print((key, " ", self.par[key]))
          # now set the check marks
          # get the menus

          menus = [i.menu() for i in self.mbar.actions()]
          items = None
          item_dict = {}
          # search for the options menu and get the contents
          for m in menus:
               if m.title() == "&Options":
                    items = m.actions()
                    for item in items:
                         if item.text() == "":
                              continue
                         else:
                              key = item.text().strip()[1:]
                              item_dict[key] = item
          # set the values for the options menu
          if item_dict != {}:
               item_dict["Convert Integer"].setChecked(self.par["convert_int"]); print(("Convert Integer = ", self.par["convert_int"]))
               item_dict["SetLimits"].setChecked(self.par["limits"]); print(("SetLimits = ", self.par["limits"]))
               item_dict["Measure"].setChecked(self.par["measure"]); print(("Measure = ", self.par["measure"]))
               item_dict["Use Limits"].setChecked(self.par["use_limits"]); print(("Use Limits = ",self.par["use_limits"]))
               item_dict["Auto Histogram Limits"].setChecked(self.par["auto_histo"]); print(("Auto Histogram Limits = ", self.par["auto_histo"]))
               item_dict['Histogram Points'].setChecked(self.par['plot_histo_points']);print(("Plot Histogram Points = ", self.par["plot_histo_points"]))
               item_dict["Filtered"].setChecked(self.par["filtered"]); print(("Filtered = ", self.par["filtered"]))
               item_dict["Plot Lines"].setChecked(self.par["draw_lines"]); print(("draw_lines = ", self.par["draw_lines"]))
          self.select_data()
          # all done


    def OnSaveParameters(self,event):
          # Create a file-open dialog in the current directory
          filetypes = '*.data'
          # use a regular file dialog
          if self.par_dir == None:
              self.par_dir = os.getcwd()

          file_dlg=QtWidgets.QFileDialog.getSaveFileName(self, 'Select parameter file to save', self.dir,
                                                         filetypes)
          if file_dlg[0] != '':
               # User has selected something, get the path, set the window's title to the path
               filename = file_dlg[0]
               # store relevant file information
               print(("will save to " , filename))
          else:
               print("so, you changed your mind, I will do nothing")
               filename = None
               return
          o = open(filename, 'w')
          # write the parameters
          for key in self.par:
               o.write( "%s  =  %r\n "%(key, self.par[key] ))
          o.close()
    def UpdateStatusBar(self, event):
          if event.inaxes:
               x, y = event.xdata, event.ydata
               self.stBar2.setText( "t= " + "%e"%(x) + "  V= " + "%e"%(y))
    def OnCloseWindow(self):

          if self.histoframe != None:
               try:
                    self.histoframe.destroy()
               except:
                    print("cannot destroy histoframe !")
          if self.tsplotframe != None:
               try:
                    self.tsplotframe.destroy()
               except:
                    print("cannot destroy tsplotframe !")
          print("Closing All")
          self.parent.quit()
          #self.destroy()
          # all done
    def OnPlot(self):
          V = None
          V = self.V
          if (self.t is None) or (V is None):
               print("Nothing to plot !")
               return
          if self.par["draw_lines"]:
               self.all_plot=self.my_plot(self.t[:V.shape[0]], V)
          else:
               self.all_plot=self.my_plot(self.t[:V.shape[0]], V, '.')
          # set the proper limits
          try:
              self.axes.set_xlim(self.t[0], self.t[:V.shape[0]][-1])
          except:
              print(f'OnPlot: cannot set limits : t-shape = {self.t.shape}, V-shape = {V.shape}')
              return
          self.axes.set_xlabel('t')
          self.axes.set_ylabel('V')
          self.figure_canvas.draw()

    def OnClear(self):
          self.figure.clf()
          self.axes=self.figure.add_subplot(111)
          self.RS = RectangleSelector(self.axes, self.LineSelectCallback,
                                      drawtype='box', useblit=True,
                                      button=[1,3], # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels')
          self.figure_canvas.draw()
    
    def OnSetScanDir(self):

         print("Set data directory")
         dir_dlg=QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a directory')
         if dir_dlg != '':
             # User has selected something, get the path, set the window's title to the path
             # store relevant file information
             self.datadir=dir_dlg
         else:
             print("so, you changed your mind, I will do nothing")
             return
         print(('Data directory: ', self.datadir))
         #self.proclistload()
          
    def OnSelectChannel(self):
          pkeys = ['chan_num']
          data = {pkeys[0] : "%d"%( int( self.par["chan_num"] ) ) }
          pdlg = NumberDialog(data, self, title="Detector channel", labels = pkeys, keys = pkeys, about_txt = "Select the data channel")
          pdlg.exec_()
          # now set the new parameters
          self.par['Detector channel']=int(pdlg.data[pkeys[0]])
          pdlg.destroy()
    
    def get_time_window(self, tmin, tmax):
          # find range of time values between tmin and tmax
          nmin = int( (tmin - self.t0)/self.dt )
          nmax = min( (int( (tmax - self.t0)/self.dt ) + 1), self.nall -1 )
          return slice(nmin, nmax + 1)
    
    def select_data(self):
          # select data range to work on
          # debugging
          if self.ndata == 0:
               return
          if (self.par["tmin"] < self.tall[0]):
               self.par["tmin"] = self.tall[0]
          if (self.par["tmax"] > self.tall[-1:][0]):
               self.par["tmax"] = self.tall[-1:][0]
          print("Get Window")
          sl = self.get_time_window( self.par["tmin"], self.par["tmax"])
          # the final value
          print("Recalculate")
          self.ndata = sl.stop - sl.start
          
          self.V = self.scale[0] + self.scale[1]*self.ydata[sl]
          self.t = self.tall[sl]
          
          print("Finished recalculation")
          
    def LineSelectCallback(self, evnt_click, evnt_release):
          'evnt_click and evnt_release are the press and release events'
          start_pos = np.array([evnt_click.xdata, evnt_click.ydata])
          end_pos = np.array([evnt_release.xdata, evnt_release.ydata] )
          self.par["xmin"] = start_pos[0]
          self.par["xmax"] = end_pos[0]
          if self.par["xmin"] > self.par["xmax"]:
               xmax = self.par["xmin"]
               self.par["xmin"] = self.par["xmax"]
               self.par["xmax"] = xmax
          self.ymin = start_pos[1]
          self.ymax = end_pos[1]
          if self.ymin > self.ymax:
               ymax = self.ymin
               self.ymin = self.ymax
               self.ymax = ymax
          dpos = end_pos - start_pos
          print(("x positions    : %5.2e --> %5.2e" % (start_pos[0],  end_pos[0] )))
          print(("y positions    : %5.2e --> %5.2e" % (start_pos[1],  end_pos[1] )))
          print(("displacement :  delta-x = %5.2e   delta-y  = %5.2e " % tuple(dpos)))
          if self.par["measure"]:
               # draw
               xpath = [self.par["xmin"], self.par["xmax"], self.par["xmax"], self.par["xmin"]]
               ypath = [self.ymin, self.ymin, self.ymax, self.ymax]
               self.axes.fill(xpath, ypath, 'b', alpha = 0.5, edgecolor = 'k')
               self.axes.text(self.par["xmin"], self.ymin, "%5.2e, %5.2e"%(start_pos[0], start_pos[1]),
                              ha = 'right', va = 'top')
               self.axes.text(self.par["xmax"], self.ymax, "%5.2e, %5.2e"%(end_pos[0], end_pos[1]) )
          if self.par["limits"]:
               ymin,ymax = self.axes.get_ylim()
               xpath = [self.par["xmin"], self.par["xmax"], self.par["xmax"], self.par["xmin"]]
               ypath = [ymin, ymin, ymax, ymax]
               self.axes.fill(xpath, ypath, 'g', alpha = 0.5, edgecolor = 'g')
               # self.axes.axvspan(self.par["xmin"], self.par["xmax"], color = 'g', alpha = 0.5)
          self.figure_canvas.draw()
          
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    frame = Main_PlotFrame(app)
    frame.show()
    sys.exit(app.exec_())