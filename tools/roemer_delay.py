from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os






##Setup plotting environment
plt.style.use('science')

h = 10
w = 10

fig = plt.figure(figsize=(w,h))
ax1 = plt.subplot2grid((1,1), (0,0))

#Load data
path = '/Users/tomkimpson/Data/Tycho/RoemerDelay/'



a = 0.6 #spin parameter

#Observer Location
ObsTheta = np.pi/2
ObsPhi = 0.0

Ox = np.sin(ObsTheta)*np.cos(ObsPhi)
Oy = np.sin(ObsTheta)*np.sin(ObsPhi)
Oz = np.cos(ObsTheta)


#some constants
Msolar = 1.989e30
c = 3e8
G=6.67e-11
MBH=4.31e6*Msolar
convert_m = c**2/(G*MBH)
convert_s = convert_m * c
convert_year = 365*24*3600




def plot(f):
    data = np.loadtxt(f)

    r = data[:,16]
    theta = data[:,17]
    phi = data[:,18]
    tau = data[:,19] #this is in seconds from rk.f of orbital dynamics
    tau = tau / convert_year    

    #Convert to cartesian
    m = np.sqrt(r**2 + a**2)
    x = m*np.sin(theta)*np.cos(phi) / convert_m
    y = m*np.sin(theta)*np.sin(phi) /convert_m
    z = r*np.cos(theta) /convert_m


    #get the roemer delay
    roemer = (x*Ox + y*Oy + z*Oz) / c
    roemer = roemer / (60*60) #hours

    return tau,roemer*1e6


eccs = ['e07/','e08/','e09/']
for e in eccs:
    f1 = path+e+'kerr.txt'
    f2 = path+e+'mpd.txt'

    t1,r1 = plot(f1)
    t2,r2 = plot(f2)

    dr = r2-r1
    ax1.plot(t1,dr)


fs = 20


ax1.locator_params(axis='both', nbins=5) #set number of xticks
ax1.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers

ax1.set_ylabel(r'$\delta_{\gamma} (\Delta_{\rm R})$ [$\mu$s]',fontsize=fs)
ax1.set_xlabel(r'$\tau$ [yr]',fontsize=fs)


savepath = '/Users/tomkimpson/Dropbox/MSSL/Papers/Submitted/Paper3_A&A/figures/'
plt.savefig(savepath+'RoemerDelay.png', dpi = 300,bbox='tight')

plt.show()
