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

num_args = len(sys.argv)

if num_args ==1:
    sys.exit('Specify how many dimensions you want to plot in. e.g. python plot_orbit_and_rays 2')



d = sys.argv[1]

print ('Plotting in ',d,'dimensions')
d = int(d)


if d == 3:
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')

if d == 2:
    print ('sdsdsdsd')
    fig, ax1,= plt.subplots(1, 1, figsize=(10,10))




#Load data
path = '/Users/tomkimpson/Data/ThesisData/'

rays = glob.glob(path + 'RT/*.txt')
orbit = path + 'MPD/trajectory.txt'
targets = path + 'MPD/targets.txt'


def Format2D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel(r'$x [r_g]$',fontsize=fs)
    ax.set_ylabel(r'$y [r_g]$',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,c='r')

    #axes limits
    sq = 1200
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)




def Format3D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel(r'$x [r_g]$',fontsize=fs)
    ax.set_ylabel(r'$y [r_g]$',fontsize=fs)
    ax.set_zlabel(r'$z [r_g]$',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,0,c='r')

    #axes limits
    sq = 1200
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)
    ax.set_zlim(-sq,sq)








def plot_ray(f,plot_type):
    data = np.loadtxt(f)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    if (plot_type == 'plot'):

        if d==2:
            ax1.plot(x,y)
        if d==3:
            ax1.plot(x,y,z)


    if (plot_type == 'scatter'):
        if d ==2:
            ax1.scatter(x,y)
        if d==3:
            ax1.scatter(x,y,z)

def plot_orbit(f,plot_type):

    data = np.loadtxt(f)

    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    if (plot_type == 'plot'):

        if d==2:
            ax1.plot(x,y,c='C0')
        if d==3:
            ax1.plot(x,y,z,c='C0')


    if (plot_type == 'scatter'):
        if d==2:
            ax1.scatter(x,y)
        if d==3:
            ax1.scatter(x,y,z)


for f in rays:
    plot_ray(f,'plot')

plot_orbit(orbit,'plot')
plot_orbit(targets,'scatter')



plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20


if d==2:
    Format2D(ax1)
if d==3:
    Format3D(ax1)

#print ('jere')
plt.show()



























