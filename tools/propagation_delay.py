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
w = 20

fig = plt.figure(figsize=(w,h))

ax0 = plt.subplot2grid((3,2), (0,0), rowspan = 3)
ax1 = plt.subplot2grid((3,2), (0,1))
ax2 = plt.subplot2grid((3,2), (1,1),sharex=ax1)
ax3 = plt.subplot2grid((3,2), (2,1),sharex=ax1)



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

def plot_traj(f,ID,col):
    data = np.loadtxt(f)

    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]


    #Plot it
    if (ID == 0):
        ax0.plot(x,y,c=col)

    if ID == 1:
        ax0.scatter(x[0],y[0],c=scat_colors[0],marker=ms)
        ax0.scatter(x[-1],y[-1],c=scat_colors[1],marker=ms)
        #mid = int(len(x)/2)

        tau = data[:,-1] #proper time in seconds
        tau = tau / convert_year
        tau = tau - tau[0]


        mid = min(range(len(tau)), key=lambda i: abs(tau[i]-0.05))
        ax0.scatter(x[mid], y[mid], c=scat_colors[2],marker=ms)





        ii = np.argwhere(x < 0)
        jj = np.argmin(np.abs(y[ii]))
        ax0.scatter(x[jj],y[jj],c=scat_colors[3],marker=ms)
 

        return jj


#Load trajecotry data
#files = glob.glob(trajpath + 'trajectory_*.txt')
#targets = glob.glob(trajpath + 'targets_*.txt')

#for f in files:
 #   plot_traj(f,0)

#for t in targets:
 #   idx = plot_traj(t,1)





def load_and_plot(f,ID,c):

    data = np.load(f)
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]

    if ID == 'alpha':
        ax1.plot(t,x/hr,c=c)
        ax1.plot(t,y/hr,c=c)
        ax2.plot(t,y-x,c=c)
        print ('Length alpha =', len(t))

        #some scatter points for reference with orbital dynamics figure
        ax1.scatter(t[0],x[0]/hr,c=scat_colors[0])
        ax1.scatter(t[-1],x[-1]/hr,c=scat_colors[1])
       
        mid = int(len(x)/2)
        mid = min(range(len(t)), key=lambda i: abs(t[i]-0.05))
        print ('mid time = ', t[mid])


        ax1.scatter(t[mid],x[mid]/hr,c=scat_colors[2])
        ax1.scatter(t[idx],x[idx]/hr,c=scat_colors[3])


        #lensing
        tL = data[:,3]
        ax2.plot(t,y-tL,c=c,linestyle='--')


       # ax1.axvline(t[idx],linestyle=ls,c=gr)
       # ax2.axvline(t[idx],linestyle=ls,c=gr)
#        ax3.axvline(t[0],linestyle=ls,c=gr)
#        ax3.axvline(t[-1],linestyle=ls,c=gr)
#

    if ID == 'beta':
        ax1.plot(t,y/hr,c=c)
        ax3.plot(t,y-x,c=c)
        print ('Length beta =', len(t))










trajpath = '/Users/tomkimpson/Data/Tycho/trajectories/'
path = 'subdata/'
eccentricities = ['e09','e08', 'e07'] #, 'e07/']
eccentricities = ['e07','e08', 'e09'] #, 'e07/']

counter = 0
for e in eccentricities:

    c= 'C'+str(counter)

    #Plot the trajectory
    plot_traj(trajpath+'trajectory_'+e+'.txt',0,c)
    idx = plot_traj(trajpath+'targets_'+e+'.txt',1,c)


    print (e,idx)


    falpha = path+'alpha_e='+e+'.npy'
    fbeta = path+'beta_e='+e+'.npy'
    


    load_and_plot(falpha,'alpha',c)
    load_and_plot(fbeta,'beta',c)
    counter = counter + 1










#Plot formatting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20



#ax2.set_yscale('log')

plt.subplots_adjust(hspace=0.01)

all_axes = plt.gcf().get_axes()
for ax in all_axes:
    ax.locator_params(axis='both', nbins=5) #set number of xticks
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
 




ax0.set_xlabel(r'$x [r_g]$',fontsize=fs)
ax0.set_ylabel(r'$y [r_g]$',fontsize=fs)

ax1.set_ylabel(r'$\Delta_{\rm prop}$ [hr]',fontsize=fs)
ax2.set_ylabel(r'$ \delta_{\alpha} (\Delta_{\rm prop})$ [s]',fontsize=fs)
ax3.set_ylabel(r'$ \delta_{\beta} (\Delta_{\rm prop})$ [s]',fontsize=fs)


ax3.set_xlabel(r'$\tau$ [yr]',fontsize=fs)
   
ax0.scatter(0,0,c='k',marker='x')
savepath = '/Users/tomkimpson/Dropbox/MSSL/Papers/Submitted/Paper3_A&A/figures/'
plt.savefig(savepath+'PropDelay.png', dpi = 300,bbox='tight')
plt.show()




















