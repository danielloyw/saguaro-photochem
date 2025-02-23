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

srun = raw_input('Enter run id > ')
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


plt.rcParams['xtick.labelsize'] = 8 
plt.rcParams['ytick.labelsize'] = 8 

clrs = ['red','blue','green','magenta','cyan']
ymin = 0
ymax = 1500

npg=0
l=0
with PdfPages('moles.pdf') as pdf:
    
    while (l < nmol):
        
        npg = npg + 1
        print npg
        fig = plt.figure()  # create a figure object
        
        # upper left
        ax1 = fig.add_subplot(2, 2, 1)
        l1 = l 
        l2 = max([l1,min([l + 5,nmol-1])])
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1])))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1+1:l2+1])))+1)
        xmin = xmax/1.E8
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)
        k = 0
        while ((k<5) & (l<nmol-1)):
            l = l + 1
            ax1.semilogx(xmol[:,l],alt,color=clrs[k])
            ax1.text(xmax/20,ymax-100-k*100,name[l],color=clrs[k],fontsize='8')
            k = k + 1
            
        ax1.set_ylabel(r'Altitude (Km)',fontsize='8')
        if (l == nmol-1): 
            pdf.savefig()
            plt.close()            
            break
    
        # upper right
        ax2 = fig.add_subplot(2, 2, 2)
        l1 = l 
        l2 = max([l1,min([l + 5,nmol-1])])
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1])))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1+1:l2+1])))+1)
        xmin = xmax/1.E8
        ax2.set_xlim(xmin,xmax)
        ax2.set_ylim(ymin,ymax)
        k = 0
        while ((k<5) & (l<nmol-1)):
            l = l + 1
            ax2.semilogx(xmol[:,l],alt,color=clrs[k])
            ax2.text(xmax/20,ymax-100-k*100,name[l],color=clrs[k],fontsize='8')
            k = k + 1
            
        if (l == nmol-1): 
            ax1.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')                 
            ax2.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')            
            pdf.savefig()
            plt.close()            
            break

    
        # lower left
        ax3 = fig.add_subplot(2, 2, 3)
        l1 = l 
        l2 = max([l1,min([l + 5,nmol-1])])
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1])))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1+1:l2+1])))+1)
        xmin = xmax/1.E8
        ax3.set_xlim(xmin,xmax)
        ax3.set_ylim(ymin,ymax)
        k = 0
        while ((k<=4) & (l<nmol-1)):
            l = l + 1
            ax3.semilogx(xmol[:,l],alt,color=clrs[k])
            ax3.text(xmax/20,ymax-100-k*100,name[l],color=clrs[k],fontsize='8')
            k = k + 1
            
        ax3.set_ylabel(r'Altitude (Km)',fontsize='8')
        ax3.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')
        if (l == nmol-1):
            ax2.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')            
            pdf.savefig()
            plt.close() 
            break
        
        # lower right
        ax4 = fig.add_subplot(2, 2, 4)
        l1 = l + 1
        l2 = max([l1,min([l + 5,nmol-1])])
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1])))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(xmol[:,l1+1:l2+1])))+1)
        xmin = xmax/1.E8
        ax4.set_xlim(xmin,xmax)
        ax4.set_ylim(ymin,ymax)
        k = 0
        while ((k<=4) & (l<nmol-1)):
            l = l + 1
            ax4.semilogx(xmol[:,l],alt,color=clrs[k])
            ax4.text(xmax/20,ymax-100-k*100,name[l],color=clrs[k],fontsize='8')
            k = k + 1

        ax4.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')
        pdf.savefig()
        plt.close()
        if (l == nmol-1):
            pdf.savefig()
            plt.close()
            break