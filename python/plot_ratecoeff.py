# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:28:00 2015

@author: rogeryelle
"""
import os
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from read_rates import read_rates


eps = 1.E-50

srun = raw_input(' Enter run id > ')
srun.strip()

#sflag = raw_input('Enter 0=screen, 1=pdf > ')
#pflag = np.int(sflag)
pflag = 1
splotfile = '../runs/'+srun+'/plots/ratecoeff.pdf'

rates = read_rates('../runs/'+srun+'/output/ratecoeff.out')
alt = 1.E-5*rates["alt"]
nalt = len(alt)

title = crates["title"] 
nrct = len(title)

rct = np.asarray(rates["rct"], dtype=np.float)

plt.rcParams['xtick.labelsize'] = 8 
plt.rcParams['ytick.labelsize'] = 8 

clrs = ['red','blue','green','magenta','cyan']
ymin = 0
ymax = 300
ydel = 20
xdel = 100000

if(pflag == 1):       
    pdf = PdfPages(splotfile)
    
 
npg=0
l=0    
while (l < nrct-1):
    
    fig = plt.figure()    
    fig.suptitle(' Run ID: '+srun)

    npg = npg + 1
    print npg
    
    # upper left
    
    print 'Upper Left'
    l1 = l 
    l2 = max([l1,min([l + 5,nrct-1])])
    if (l1 == l2):
        xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1]+eps)))+1)
    else:    
        xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1:l2+1]+eps)))+1)
    xmin = xmax/1.E8
    
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    k = 0
    while ((k<5) & (l<nrct-1)):
        ax1.semilogx(rct[:,l],alt,color=clrs[k])
        ax1.text(xmax/xdel,ymax-ydel-k*ydel," ".join(title[l].split()),color=clrs[k],fontsize='8')
        l = l + 1
        k = k + 1
#        print title[l]
        
    ax1.set_ylabel(r'Altitude (Km)',fontsize='8')
    if (l == nrct-1): 
        ax1.set_xlabel(r'Rate Coefficient (cm$^{3}$s$^{-1}$)',fontsize='8')                          


    # upper right
    print 'Upper Right'
    if(l<nrct-1):
        l1 = l 
        l2 = max([l1,min([l + 5,nrct-1])])   
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1]+eps)))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1:l2+1]+eps)))+1)
        xmin = xmax/1.E8
        
        ax2 = fig.add_subplot(2, 2, 2)
        ax2.set_xlim(xmin,xmax)
        ax2.set_ylim(ymin,ymax)
        k = 0
        while ((k<5) & (l<nrct-1)):
            ax2.semilogx(rct[:,l],alt,color=clrs[k])
            ax2.text(xmax/xdel,ymax-ydel-k*ydel," ".join(title[l].split()),color=clrs[k],fontsize='8')
            k = k + 1
            l = l + 1
#            print title[l]
            
            if (l == nrct-1): 
                ax1.set_xlabel(r'Rate Coefficient (cm$^{3}$s$^{-1}$)',fontsize='8')                 
                ax2.set_xlabel(r'Rate Coefficient (cm$^{3}$s$^{-1}$)',fontsize='8')                      
    

    # lower left
    print 'Lower Left'
    if(l<nrct-1):
        l1 = l 
        l2 = max([l1,min([l + 5,nrct-1])])
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1]+eps)))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1:l2+1]+eps)))+1)
        xmin = xmax/1.E8
        
        ax3 = fig.add_subplot(2, 2, 3)
        ax3.set_xlim(xmin,xmax)
        ax3.set_ylim(ymin,ymax)
        
        k = 0
        while ((k<=4) & (l<nrct-1)):
            ax3.semilogx(rct[:,l],alt,color=clrs[k])
            ax3.text(xmax/xdel,ymax-ydel-k*ydel," ".join(title[l].split()),color=clrs[k],fontsize='8')
            k = k + 1
            l = l + 1
#            print title[l]
            
        ax3.set_ylabel(r'Altitude (Km)',fontsize='8')
        ax3.set_xlabel(r'Rate Coefficient (cm$^{3}$s$^{-1}$)',fontsize='8')
    
        if (l == nrct-1):
            ax2.set_xlabel(r'Rate Coefficient (cm$^{3}$s$^{-1}$)',fontsize='8')            
    
    # lower right
    print 'Lower Right'
    if(l<nrct-1):
        l1 = l
        l2 = max([l1,min([l + 5,nrct-1])])    
        if (l1 == l2):
            xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1]+eps)))+1)
        else:    
            xmax = math.pow(10,int(math.log10(np.amax(rct[:,l1:l2+1]+eps)))+1)
        xmin = xmax/1.E8
        
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.set_xlim(xmin,xmax)
        ax4.set_ylim(ymin,ymax)
        
        k = 0
        while ((k<=4) & (l<nrct-1)):
            ax4.semilogx(rct[:,l],alt,color=clrs[k])
            ax4.text(xmax/xdel,ymax-ydel-k*ydel," ".join(title[l].split()),color=clrs[k],fontsize='8')
            k = k + 1
            l = l + 1
#            print title[l]
    
        ax4.set_xlabel(r'Rate Coefficient (cm$^{3}$s$^{-1}$)',fontsize='8')
    
    #
    #  Close this page
    #
    
    if(pflag == 1):
        pdf.savefig(fig)
        plt.close()
    else:
        plt.show()  
    
#
# close the pdf file for this molecule
#
    
if(pflag == 1):           
    pdf.close()              
  