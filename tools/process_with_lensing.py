from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os
from scipy import interpolate


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


def calculate_tPK(x1,x2):


        dr = x2-x1
        mag = np.linalg.norm(dr)
        u = dr / mag
        dp = np.dot(x1,u)
        tPK = mag + 2*np.log((mag+dp+np.linalg.norm(x2))/(dp+np.linalg.norm(x1)))

        return tPK

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
        tPK = calculate_tPK(x1,x2)

        #lensing
        #interpolate - this will be used for lensing corrections
        
        if min(xH) < 0:
            fxy = interpolate.interp1d(xH, yH)
            fxz = interpolate.interp1d(xH, zH)


            xc = [0,fxy(0),fxz(0)] #coordinates of lensing plane point
        

            tL1 = calculate_tPK(x1,xc)
            tL2 = calculate_tPK(xc,x2)
            tL = tL1 + tL2
 #           print ('lensing')
        else:
            tL = tPK
  #          print ('nolensing')
    

 #       print (tPK/convert_s, tL/convert_s)
        return tau/(convert_year), tf/convert_s,tPK/convert_s, tL/convert_s

    else:
        return tau/(convert_year),tf/convert_s, PK,0



def process(e,typ,PK):

    rays = glob.glob(path+typ+e+'*.txt')

    #Create some arrays to hold the t/tau data
    xx = []
    yy = []


    zz = [] #an extra array to hold tPK, as calculated from Schwarzchild soln
    aa = [] #an extra array to hold tPK_lensing, as calculated from Schwarzchild soln

    
    for f in rays:
        tau,t,tPK,tL = extract(f,PK)       
        xx.extend([tau])
        yy.extend([t])
        zz.extend([tPK])
        aa.extend([tL])

    #sort it and shift it 
    Z = [x for _,x in sorted(zip(xx,yy))]
    yy  = Z
    yy = yy-yy[0]
    

    if PK == 1:
        Z = [x for _,x in sorted(zip(xx,zz))]
        zz  = Z
  #      print ('TPK',zz[0])
   #     print (zz)
        zz =zz - zz[0]

        Z = [x for _,x in sorted(zip(xx,aa))]
        aa  = Z
    #    print ('TL',aa[0])
     #  print (aa)
        aa =aa - aa[0]


    xx = sorted(xx) 


    return np.array(xx),np.array(yy),np.array(zz),np.array(aa) #zz/aa are only meaningful if PK=1



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
    print (estr)
  
    ls = linestyles[counter]

    #Do stuff with schwarzchild and PK, which deffo have the same tau
    tau1,tSCH,tPK1,tL1 = process(e,types[0],1)
    tau1 = tau1-tau1[0]

    f1 = 'subdata/alpha_e='+estr+'.npy'
    
    s1 = np.zeros((len(tau1),4))
    s1[:,0] = tau1
    s1[:,1] = tPK1
    s1[:,2] = tSCH
    s1[:,3] = tL1

    np.save(f1, s1)


    #Now compare the schwarzchild and kerr solutions
    #these might now have the same tau
    tau2,tKER,tPK2,tL2 = process(e,types[1],0)
    tau2 = tau2-tau2[0]

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


    counter = counter + 1




print ('Completed')
