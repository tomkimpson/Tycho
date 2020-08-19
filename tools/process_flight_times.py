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
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, figsize=(10,10),sharex=True)





#Load data
path = '/Users/tomkimpson/Data/Tycho/PropDelay/'


#some constants
Msolar = 1.989e30
c = 3e8
G=6.67e-11
MBH=4.31e6*Msolar
convert_m = c**2/(G*MBH)
convert_s = convert_m * c
convert_year = 365*24*3600


hr = 60*60


def extract(f,PK):

    data = np.loadtxt(f)

    t=data[:,5]
    tf = t[-1] - t[0] #flight time
    tau = data[0,9]
    


    if PK ==1:
        #Do some calculations to get the post-Keplerian flight time
        

        rBL = data[:,10]

        theta = data[:,11]
        phi = data[:,12]
        rh = rBL - 1

        xH = rh*np.sin(theta)*np.cos(phi)
        yH = rh*np.sin(theta)*np.sin(phi)
        zH = rh*np.cos(theta)


        x1 = np.array([xH[-1], yH[-1], zH[-1]])
        x2 = np.array([xH[0], yH[0], zH[0]])


        dr = x2-x1
        mag = np.linalg.norm(dr)
        u = dr / mag
        
        dp = np.dot(x1,u)


        tPK = mag + 2*np.log((mag+dp+np.linalg.norm(x2))/(dp+np.linalg.norm(x1)))

        return tau/(convert_year), tf/convert_s,tPK/convert_s

    else:
        return tau/(convert_year),tf/convert_s, PK



def process(e,typ,PK):

    rays = glob.glob(path+typ+e+'*.txt')

    #Create some arrays to hold the t/tau data
    xx = []
    yy = []


    zz = [] #an extra array to hold tPK, as calculated from Schwarzchild soln

    
    for f in rays:
        tau,t,tPK = extract(f,PK)       
        xx.extend([tau])
        yy.extend([t])
        zz.extend([tPK])


    #sort it and shift it 
    Z = [x for _,x in sorted(zip(xx,yy))]
    yy  = Z
    yy = yy-yy[0]
    

    if PK == 1:
        Z = [x for _,x in sorted(zip(xx,zz))]
        zz  = Z
        zz =zz - zz[0]

    xx = sorted(xx) 
    return np.array(xx),np.array(yy),np.array(zz) #zz is only meaningful if PK=1



#eccentricities = ['e07/', 'e08/', 'e09/']


types = ['schwarz/', 'kerr/']


eccentricities = ['e09/','e08/', 'e07/'] #, 'e07/']
linestyles = ['solid','dashed', 'dotted']
counter = 0

#eccentricities = ['e07/'] #, 'e07/']

counter = 0
#for e in eccentricities:




for e in eccentricities:


    estr= str(e[:-1])

  
    ls = linestyles[counter]

    #Do stuff with schwarzchild and PK, which deffo have the same tau
    tau1,tSCH,tPK1 = process(e,types[0],1)
    tau1 = tau1-tau1[0]
    ax1.plot(tau1,tPK1/hr,linestyle=ls)
    ax1.plot(tau1,tSCH/hr,linestyle=ls)

    f1 = 'subdata/alpha_e='+estr+'.npy'
 
    s1 = np.zeros((len(tau1),3))
    s1[:,0] = tau1
    s1[:,1] = tPK1
    s1[:,2] = tSCH

    np.save(f1, s1)


    alpha = tSCH - tPK1
    ax2.plot(tau1,np.abs(alpha),c='C0',linestyle=ls)
#    ax2.scatter(tau1,np.abs(alpha),c='C0',linestyle=ls)
 #   ax2.scatter(tau1,alpha)





    #Now compare the schwarzchild and kerr solutions
    #these might now have the same tau
    tau2,tKER,tPK2 = process(e,types[1],0)
    tau2 = tau2-tau2[0]
    ax1.plot(tau2,tKER/hr,linestyle=ls)
#    ax1.scatter(tau2,tKER,c='r')



    tau_crop = []
    tSCH_crop= []
    tKER_crop= []
    for i in range(len(tau1)):
        if tau1[i] in tau2:
            j = np.argwhere(tau1[i] == tau2)
            j = j[0][0]
            tau_crop.extend([tau1[i]])
            tSCH_crop.extend([tSCH[i]])
            tKER_crop.extend([tKER[j]])
            
 #           print (i,j, tau1[i],tau2[j])
#
    
    tau_crop = np.array(tau_crop)
    tSCH_crop = np.array(tSCH_crop)
    tKER_crop = np.array(tKER_crop)



    f2 = 'subdata/beta_e='+estr+'.npy'
 
    s2 = np.zeros((len(tau_crop),3))
    s2[:,0] = tau_crop
    s2[:,1] = tSCH_crop
    s2[:,2] = tKER_crop

    np.save(f2, s2)


    beta = tKER_crop - tSCH_crop
    ax3.plot(tau_crop, np.abs(beta),c='C1',linestyle=ls)



    print ('Minima:', min(np.abs(alpha[1:])), min(np.abs(beta[1:])))

 #   check = np.array_equal(tau1,tau2)
 #   if check == False:
 #       sys.exit('Exited as the two tau arrays are different')





    counter = counter + 1



plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20



#ax2.set_yscale('log')

plt.subplots_adjust(hspace=0.01)

all_axes = plt.gcf().get_axes()
for ax in all_axes:
    ax.locator_params(axis='both', nbins=5) #set number of xticks
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
 


ax1.set_ylabel(r'$\Delta_{\rm prop}$ [hr]',fontsize=fs)
ax2.set_ylabel(r'$ \delta_{\alpha} (\Delta_{\rm prop})$ [s]',fontsize=fs)
ax3.set_ylabel(r'$ \delta_{\beta} (\Delta_{\rm prop})$ [s]',fontsize=fs)


ax3.set_xlabel(r'$\tau$ [yr]',fontsize=fs)
   

savepath = '/Users/tomkimpson/Dropbox/MSSL/Papers/Submitted/Paper3_A&A/figures/'
plt.savefig(savepath+'PropDelay.png', dpi = 300,bbox='tight')
plt.show()




















