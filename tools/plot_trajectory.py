from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os





d = 3


#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig = plt.figure(figsize=(20,10))
    ax1 = plt.subplot2grid((1,2), (0,0))
    ax2 = plt.subplot2grid((1,2), (0,1))



#Load data
path = '/Users/tomkimpson/Data/Tycho/'
files = glob.glob(path + 'trajectory_*.txt')



def plot(f):
    data = np.loadtxt(f)

    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]


    r = x**2 + y**2 + z**2
    idx = np.argmin(r)
    print (idx)



    #Plot it


    if (d == 3):
        ax1.plot(x,y,z)
        ax1.scatter(x[0],y[0],z[0], c='g')
        ax1.scatter(x[-1],y[-1],z[-1], c='r')
        ax1.scatter(0,0,0, c='k')

 
        perix = [0,x[idx]]
        periy = [0,y[idx]]
        periz = [0,z[idx]]
        ax1.plot(perix,periy,periz,c='C1')

        ax1.scatter(x[idx],y[idx],z[idx], c='C1')


        limit = max(max(x),max(y),max(z))
        ax1.set_xlim(-limit,+limit)
        ax1.set_ylim(-limit,+limit)
        ax1.set_zlim(-limit,+limit)

    if (d == 2):
        ax1.plot(x,y)
        ax2.plot(x,z)


for f in files:
    plot(f)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

plt.show()
