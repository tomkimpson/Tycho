from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os
from scipy.signal import argrelextrema
from scipy import interpolate


#Setup plotting environment
plt.style.use('science')
fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True, figsize=(10,10))



#Load data
path = '/Users/tomkimpson/Data/Tycho/'
PlotFile = path + 'PK.txt'

hr = 3600
yr = 3600*24*365


def get_data(f):

    data = np.loadtxt(f,skiprows=1)
    t = data[:,0] / yr
    roemer = data[:,2] #/ hr

    return t,roemer




def process(ls):


    #Define the files
    F2 = path+'TimeDelay_PK.txt'

    #get the data
    x2,y2 = get_data(F2)


    #plot the roemer delay
    ax1.plot(x2,y2/hr)



#Format
fs = 20

linestyles = ['solid','dashed', 'dotted']
counter = 0

ls = linestyles[counter]
process(ls)

all_axes = plt.gcf().get_axes()
for ax in all_axes:
    ax.locator_params(axis='both', nbins=5) #set number of xticks
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
    

plt.setp(ax1.get_xticklabels(),visible=False)


ax2.set_yscale('log')
plt.subplots_adjust(hspace = 0.01)
#Label axes

fontsize=fs

ax1.set_ylabel(r'$\Delta_{\rm E}$ [hr]',fontsize=fs)
ax2.set_ylabel(r'$|\delta \Delta_{\rm E}|$ [s]',fontsize=fs)
ax2.set_xlabel(r'$\tau$ [yr]',fontsize=fs)

ax2.set_ylim(1e-7)

plt.show()


