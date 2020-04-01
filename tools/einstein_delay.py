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


a = 7606
e = 0.88

def get_data(f):

    data = np.loadtxt(f)
    t = data[:,0]
    E1 = data[:,1]
    E2 = data[:,2]
    r = data[:,3]


    ax1.plot(t,E1)
    ax1.plot(t,E2,c='r')



    EA = np.arccos((1-r/a)/e)

    LHS = EA - e*np.sin(EA)
    RHS = 2*np.pi*t

    ax2.plot(t,np.sin(LHS) % 2*np.pi)
    ax2.plot(t,np.sin(RHS))

#Format
fs = 20

get_data(PlotFile)

all_axes = plt.gcf().get_axes()
for ax in all_axes:
    ax.locator_params(axis='both', nbins=5) #set number of xticks
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
    

#Label axes

fontsize=fs

ax1.set_ylabel(r'$\Delta_{\rm E}$ [hr]',fontsize=fs)



plt.show()

