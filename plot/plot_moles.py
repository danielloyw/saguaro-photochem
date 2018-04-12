# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 17:51:20 2014

@author: rogeryelle
"""

from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from read_nmolecules import read_nmolecules
from read_atm import read_atm

srun = input('Enter run id > ')
srun.strip()  

satmfile = '../runs/'+srun+'/output/atm1D.out'    
atm = read_atm(satmfile)
den = atm["den"]
alt = 1.E-5*atm["alt"]
name = atm["name"]
map(str.strip, name)

nmol = len(name)
nalt = len(alt)
xmol = np.zeros((nalt,nmol))

for m in range(1,nmol):
   xmol[:,m] = den[:,m]/den[:,0]


clrs = ['red','blue','green','magenta','cyan']
ymin = 0
ymax = 250
xmin = 1.E4

l=0
with PdfPages('../runs/'+srun+'/plots/moles.pdf') as pdf:
    
    while (l <= nmol-1):
        fig = plt.figure()  # create a figure object
        
        for subplot_i in range(1,5):
            ax = fig.add_subplot(2, 2, subplot_i)
            l1 = l
            l2 = max([l1,min([l + 5,nmol-1])])
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1:l2+1])))+1)
            xmin = min(xmin,xmax/1.E8)
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)
            k = 0
            while ((k<5) and (l<=nmol-1)):
                ax.semilogx(xmol[:,l],alt,color=clrs[k])
                ax.text(xmax/100,ymax-20-k*20,name[l],color=clrs[k],fontsize='8')
                l = l + 1
                k = k + 1
            if (subplot_i == 1) or (subplot_i == 3):
                ax.set_ylabel(r'Altitude (Km)')
            if (subplot_i == 3) or (subplot_i == 4):
                ax.set_xlabel(r'Mole Fraction (V/V)')
            if (l == nmol): 
                break
        pdf.savefig()
        plt.close()
