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
path = '/Users/tomkimpson/Data/Tycho/EinsteinDelay/'
PlotFile = path + 'PK.txt'

hr = 3600
yr = 3600*24*365


def get_data(f):

    data = np.loadtxt(f,skiprows=1)
    t = data[:,0] / yr
    Einstein = data[:,1] #/ hr

    return t,Einstein




def process(e,ls):


    #Define the files
    F1 = path+e+'TimeDelay_PK.txt'
    F2 = path+e+'TimeDelay_SCH.txt'
    F3 = path+e+'TimeDelay_KER.txt'
    F4 = path+e+'TimeDelay_MPD.txt'
    
    print (F1)


    #get the data
    x1,y1 = get_data(F1)
    x2,y2 = get_data(F2)
    x3,y3 = get_data(F3)
    x4,y4 = get_data(F4)


    #plot the Einstein delay
    ax1.plot(x1,y1/hr)
    ax1.plot(x2,y2/hr)
    ax1.plot(x3,y3/hr)
    ax1.plot(x4,y4/hr)

    #get the differences
    delta_alpha = np.abs(y2 - y1)
    delta_beta = np.abs(y3 - y2)
    delta_gamma = np.abs(y4 - y3)

    #Plot the differnces
    ax2.plot(x1,delta_alpha, c='C0',linestyle=ls)
    ax2.plot(x1,delta_beta, c='C1',linestyle=ls)
    ax2.plot(x1,delta_gamma, c='C2',linestyle=ls)

    print ('Max  alpha = ', max(delta_alpha))
    print ('Max  beta = ', max(delta_beta))
    print ('Max  gamma = ', max(delta_gamma))
    print ('-----')


#Format
fs = 20

eccentricities = ['e09/','e08/', 'e07/'] #, 'e07/']
linestyles = ['solid','dashed', 'dotted']
counter = 0

for e in eccentricities:
    ls = linestyles[counter]
    process(e,ls)
    counter = counter + 1
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

savepath = '/Users/tomkimpson/Dropbox/MSSL/Papers/Submitted/Paper3_A&A/figures/'
plt.savefig(savepath+'EinsteinDelay.png', dpi = 300,bbox='tight')
plt.show()


