# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 17:51:20 2014

@author: rogeryelle
"""

from pylab import *
from scipy.constants import G, Boltzmann as Kb, m_u
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
from matplotlib import rc
import matplotlib.cm as cm
from read_tchem1D import *


wav, name, crs_tot, crs_ion, nbrnch, reaction, brnch = read_photoA('../data/photons/photoA.dat')
nmol = len(name)
#wav = photoA["wav"]
#nwav=wav.size
#name = photoA["name"]
#nmol=len(name)
#crs_tot = photoA["crs_tot"]
#crs_ion = photoA["crs_ion"]
#brnch = photoA["brnch"]
#nbrnch = np.size(brnch,1)
#crs_brnch = zeros((nwav,nbrnch), dtype=float)
   
#uniqname = np.unique(name)

#for uname in uniqname:
#    indx = [i for i, x in enumerate(name) if x == uname]
#    crs_brnch[:,indx]=brnch[:,indx]*crs_tot[:,"uname"]

crs_brnch=np.zeros_like(brnch)
nb = 0
for nm in xrange(0,nmol):
    for ny in xrange(0,nbrnch[nm]):
        crs_brnch[:,nb]=crs_tot[:,nm]*brnch[:,nb]
        nb = nb + 1
        
       
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.subplots_adjust(hspace=0.4)

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('test2png.png',dpi=100)
        
nb = 0        
for nm in xrange(0,nmol):
    if(nm ==0):
        plt.subplot(221)
        plt.ylabel(r'Cross Section (cm$^2$)')
    elif(nm==1):
        plt.subplot(222)
    elif(nm==2):
        plt.subplot(223)
        plt.xlabel(r'Wavelength (\AA )')
        plt.ylabel(r'Cross Section (cm$^2$)')
    elif(nm==3):
        plt.subplot(224)
        plt.xlabel(r'Wavelength (\AA )')
        
    plt.grid(True)
    plt.semilogy(wav,crs_tot[:,nm])
    plt.axis([0, 800, 1.E-24, 1.E-16])    
    for nx in xrange(0,nbrnch[nm]):
        plt.semilogy(wav,crs_brnch[:,nb])
        plt.text(400,1E-19*pow(10,-0.3*nx),reaction[nb].strip(),fontsize=6)
        nb = nb + 1
    plt.text(100,1.E-17,name[nm])    
            
plt.show()

       