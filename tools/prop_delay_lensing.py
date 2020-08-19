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



hr = 60*60

#some constants
Msolar = 1.989e30
c = 3e8
G=6.67e-11
MBH=4.31e6*Msolar
convert_m = c**2/(G*MBH)
convert_s = convert_m * c
convert_year = 365*24*3600


ms = 'o'
ls = '--'
gr = '0.5'


scat_colors = ['turquoise', 'purple','#8c564b', 'r']


def load_and_plot(f,ID,c):

    data = np.load(f)
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]

    if ID == 'alpha':
        ax1.plot(t,y-x,c=c,linestyle='--')
        #get the lensing correction
        tL = data[:,3]
        ax1.plot(t,y-tL,c=c)


trajpath = '/Users/tomkimpson/Data/Tycho/trajectories/'
path = 'subdata/'
eccentricities = ['e09','e08', 'e07'] #, 'e07/']
counter = 0
for e in eccentricities:

    c= 'C'+str(counter)


    falpha = path+'alpha_e='+e+'.npy'

    load_and_plot(falpha,'alpha',c)
    counter = counter + 1






#Plot formatting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20
 

ax1.locator_params(axis='both', nbins=5) #set number of xticks
ax1.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
 
ax1.set_ylabel(r'$ \delta_{\alpha} (\Delta_{\rm prop})$ [s]',fontsize=fs)
ax1.set_xlabel(r'$\tau$ [yr]',fontsize=fs)



savepath = '/Users/tomkimpson/Dropbox/MSSL/Papers/Submitted/Paper3_A&A/figures/'
plt.savefig(savepath+'PropDelay_lensing.png', dpi = 300,bbox='tight')
plt.show()




















