# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 13:01:14 2015

@author: rogeryelle
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from read_atm import read_atm
from read_data import read_data

# plot to screen or pdf

pflag = 1

# get run id

srun = raw_input('Enter run id > ')
srun.strip()

# output file name

if(pflag == 1):
    splotfile = '../runs/'+srun+'/plots/data_plot.pdf'

# list of molecular species to be plotted

molec=['C2H2','C2H4','C2H6','CH2CCH2','CH3CCH','C3H6','C3H8','C4H2','C4H6', \
       'C5H4','C6H2','C6H6', 'C7H4','C7H8','C8H2','HCN','HNC','HC3N','CH3CN',    \
       'C3H3N','C3H5N','C4H3N','C4H5N','C5H5N','C6H3N','C6H7N','CH2NH','CH3NH2', \
       'NH3','N2H4','CO','CO2','H2O','CH3OH']

#molec = ['CO2']       
nmol = len(molec)       

# read atmosphere

atm = read_atm('../runs/'+srun+'/output/atm1D.out')

# extract atmosphere data from dictionary

den = atm["den"]
alt = 1.E-5*atm["alt"]
name = atm["name"]
map(str.strip, name)
nalt = len(alt)

# read observation files

# need to be able to pull all the data for a given molecule


# extract model mole fractions

xmol = np.zeros((nalt,nmol))
for l in xrange(0,nmol):
    xmol[:,l] = den[:,name.index(molec[l])]/den[:,0]
    
# Make plots

plt.rcParams['xtick.labelsize'] = 8 
plt.rcParams['ytick.labelsize'] = 8 

###
###  add more colors
###

clrs = ['black','red','blue','green','magenta','cyan','yellow']
ymin = 0
ymax = 1500
xrng = 1.E6

if(pflag == 1):       
    pdf = PdfPages(splotfile)
    
 
npg=0
l=0    
while (l < nmol):
    
    # set up figure
    
    fig = plt.figure()    
    fig.suptitle(' Run ID: '+srun)
    npg = npg + 1
    print npg

    #        
    # upper left
    #

    #  retrieve observations
    
    data = read_data(molec[l])
    cirs = data["cirs"]
    uvis = data["uvis"]
    inms = data["inms"]
    
    # set up panel
    
    xmx_data = max(max(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
    xmx_modl = max(xmol[:,l])
    xmax = max([xmx_data,xmx_modl])
    xmax = np.power(10,int(np.log10(xmax)))
    xmax = min(xmax,1.0)
    
    xmn_data = min(min(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
    xmn_modl = min(xmol[:,l])
    xmin = min([xmn_data,xmn_modl])
    xmin = np.power(10,int(np.log10(xmin))-1)
    xmin = max(xmin,1.E-12)
    
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_xscale('log')
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)    
    
    #  plot model
    
    ax1.plot(xmol[:,l],alt,'-k',linewidth=1) 
    
    # plot CIRS data

    ndat_cirs = cirs["ndat"]
    alt_cirs = cirs["alt"]
    xmol_cirs = cirs["xmol"]
    xmin_cirs = cirs["xmin"]
    xmax_cirs = cirs["xmax"]
    for x in xrange(0,len(ndat_cirs)):
        ax1.plot(xmol_cirs[x],alt_cirs[x],linestyle='-',marker=' ', \
        color=clrs[(x % 6)+1])  

    # plot UVIS data

    uvis = data["uvis"]
    ndat_uvis = uvis["ndat"]
    alt_uvis = uvis["alt"]
    xmol_uvis = uvis["xmol"]
    xmin_uvis = uvis["xmin"]
    xmax_uvis = uvis["xmax"]
    for x in xrange(0,len(ndat_uvis)):
        ax1.plot(xmol_uvis[x],alt_uvis[x],linestyle='--',marker=' ', \
        color=clrs[(x % 6)+1])  

    # plot INMS data    

    inms = data["inms"]
    ndat_inms = inms["ndat"]
    alt_inms = inms["alt"]
    xmol_inms = inms["xmol"]
    xmin_inms = inms["xmin"]
    xmax_inms = inms["xmax"]
    for x in xrange(0,len(ndat_inms)):
        ax1.scatter(xmol_inms[x],alt_inms[x],marker='o',color=clrs[(x % 6)+1])  

    # write axes labels
    
    ax1.text(2*xmin,ymax-100,molec[l],color=clrs[0],fontsize='8')        
    ax1.set_ylabel(r'Altitude (Km)',fontsize='8')
    l = l + 1
    if (l == nmol-1): 
        ax1.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')                          


    # upper right

    if(l<nmol-1):

        #  retrieve observations
        
        data = read_data(molec[l])
        cirs = data["cirs"]
        uvis = data["uvis"]
        inms = data["inms"]
        
        # set up panel
    
        xmx_data = max(max(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
        xmx_modl = max(xmol[:,l])
        xmax = max([xmx_data,xmx_modl])
        xmax = np.power(10,int(np.log10(xmax)))
        xmax = min(xmax,1.0)
        
        xmn_data = min(min(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
        xmn_modl = min(xmol[:,l])
        xmin = min([xmn_data,xmn_modl])
        xmin = np.power(10,int(np.log10(xmin))-1)
        xmin = max(xmin,1.E-12)
        
        ax2 = fig.add_subplot(2, 2, 2)
        ax2.set_xscale('log')
        ax2.set_xlim(xmin,xmax)
        ax2.set_ylim(ymin,ymax)    
    
        # plot model
        
        ax2.plot(xmol[:,l],alt,'-k',linewidth=1)
                
        # plot CIRS data
    
        cirs = data["cirs"]
        ndat_cirs = cirs["ndat"]
        alt_cirs = cirs["alt"]
        xmol_cirs = cirs["xmol"]
        xmin_cirs = cirs["xmin"]
        xmax_cirs = cirs["xmax"]
        for x in xrange(0,len(ndat_cirs)):
            ax2.semilogx(xmol_cirs[x],alt_cirs[x],color=clrs[(x % 6)+1])  
    
        # plot UVIS data
    
        uvis = data["uvis"]
        ndat_uvis = uvis["ndat"]
        alt_uvis = uvis["alt"]
        xmol_uvis = uvis["xmol"]
        xmin_uvis = uvis["xmin"]
        xmax_uvis = uvis["xmax"]
        for x in xrange(0,len(ndat_uvis)):
            ax2.semilogx(xmol_uvis[x],alt_uvis[x],linestyle='--',marker=' ', \
            color=clrs[(x % 6)+1])  
    
        # plot INMS data    
    
        inms = data["inms"]
        ndat_inms = inms["ndat"]
        alt_inms = inms["alt"]
        xmol_inms = inms["xmol"]
        xmin_inms = inms["xmin"]
        xmax_inms = inms["xmax"]
        
        for x in xrange(0,len(ndat_inms)):
            ax2.scatter(xmol_inms[x],alt_inms[x],linestyle=':',marker='o', \
            color=clrs[(x % 6)+1])  

        # write axes labels
            
        ax2.text(2*xmin,ymax-100,molec[l],color=clrs[0],fontsize='8')
        l = l + 1
        if (l == nmol-1): 
            ax1.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')                 
            ax2.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')                      
    

    # lower left

    if(l<nmol-1):

        #  retrieve observations
        
        data = read_data(molec[l])
        cirs = data["cirs"]
        uvis = data["uvis"]
        inms = data["inms"]
        
        # set up panel
    
        xmx_data = max(max(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
        xmx_modl = max(xmol[:,l])
        xmax = max([xmx_data,xmx_modl])
        xmax = np.power(10,int(np.log10(xmax)))
        xmax = min(xmax,1.0)
        
        xmn_data = min(min(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
        xmn_modl = min(xmol[:,l])
        xmin = min([xmn_data,xmn_modl])
        xmin = np.power(10,int(np.log10(xmin))-1)
        xmin = max(xmin,1.E-12)
         
        ax3 = fig.add_subplot(2, 2, 3)
        ax3.set_xscale('log')
        ax3.set_xlim(xmin,xmax)
        ax3.set_ylim(ymin,ymax)    
        
        #  plot model
        
        ax3.semilogx(xmol[:,l],alt,color=clrs[0])
                
        # plot CIRS data
    
        cirs = data["cirs"]
        ndat_cirs = cirs["ndat"]
        alt_cirs = cirs["alt"]
        xmol_cirs = cirs["xmol"]
        xmin_cirs = cirs["xmin"]
        xmax_cirs = cirs["xmax"]
        for x in xrange(0,len(ndat_cirs)):
            ax3.semilogx(xmol_cirs[x],alt_cirs[x],color=clrs[(x % 6)+1])  
    
        # plot UVIS data
    
        uvis = data["uvis"]
        ndat_uvis = uvis["ndat"]
        alt_uvis = uvis["alt"]
        xmol_uvis = uvis["xmol"]
        xmin_uvis = uvis["xmin"]
        xmax_uvis = uvis["xmax"]
        for x in xrange(0,len(ndat_uvis)):
            ax3.semilogx(xmol_uvis[x],alt_uvis[x],linestyle='--', \
            color=clrs[(x % 6)+1])  
    
        # plot INMS data    
    
        inms = data["inms"]
        ndat_inms = inms["ndat"]
        alt_inms = inms["alt"]
        xmol_inms = inms["xmol"]
        xmin_inms = inms["xmin"]
        xmax_inms = inms["xmax"]
        
        for x in xrange(0,len(ndat_inms)):
            ax3.scatter(xmol_inms[x],alt_inms[x],color=clrs[(x % 6)+1])  

        # write axes labels
            
        ax3.text(2*xmin,ymax-100,molec[l],color=clrs[0],fontsize='8')
        ax3.set_ylabel(r'Altitude (Km)',fontsize='8')
        ax3.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')    
        l = l + 1            
        if (l == nmol-1):
            ax2.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')            

    
    # lower right

    if(l<nmol-1):

        #  retrieve observations
        
        data = read_data(molec[l])
        cirs = data["cirs"]
        uvis = data["uvis"]
        inms = data["inms"]
        
        # set up panel
    
        xmx_data = max(max(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
        xmx_modl = max(xmol[:,l])
        xmax = max([xmx_data,xmx_modl])
        xmax = np.power(10,int(np.log10(xmax)))
        xmax = min(xmax,1.0)
        
        xmn_data = min(min(cirs["xmol"]+uvis["xmol"]+inms["xmol"]))
        xmn_modl = min(xmol[:,l])
        xmin = min([xmn_data,xmn_modl])
        xmin = np.power(10,int(np.log10(xmin))-1)
        xmin = max(xmin,1.E-12)
         
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.set_xscale('log')
        ax4.set_xlim(xmin,xmax)
        ax4.set_ylim(ymin,ymax)    
        
        #  plot model

        ax4.semilogx(xmol[:,l],alt,color=clrs[0])
        
        # plot CIRS data
    
        cirs = data["cirs"]
        ndat_cirs = cirs["ndat"]
        alt_cirs = cirs["alt"]
        xmol_cirs = cirs["xmol"]
        xmin_cirs = cirs["xmin"]
        xmax_cirs = cirs["xmax"]
        for x in xrange(0,len(ndat_cirs)):
            ax4.semilogx(xmol_cirs[x],alt_cirs[x],color=clrs[(x % 6)+1])  
    
        # plot UVIS data
    
        uvis = data["uvis"]
        ndat_uvis = uvis["ndat"]
        alt_uvis = uvis["alt"]
        xmol_uvis = uvis["xmol"]
        xmin_uvis = uvis["xmin"]
        xmax_uvis = uvis["xmax"]
        for x in xrange(0,len(ndat_uvis)):
            ax4.semilogx(xmol_uvis[x],alt_uvis[x],linestyle='--', \
            color=clrs[(x % 6)+1])  
    
        # plot INMS data    
    
        inms = data["inms"]
        ndat_inms = inms["ndat"]
        alt_inms = inms["alt"]
        xmol_inms = inms["xmol"]
        xmin_inms = inms["xmin"]
        xmax_inms = inms["xmax"]
        
        for x in xrange(0,len(ndat_inms)):
            ax4.scatter(xmol_inms[x],alt_inms[x],color=clrs[(x % 6)+1])  

        # write axes labels
            
        ax4.text(2*xmin,ymax-100,molec[l],color=clrs[0],fontsize='8')
        ax4.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')
        l = l + 1

    
    #  Close this page
    
    if(pflag == 1):
        pdf.savefig(fig)
        plt.close()
    else:
        plt.show()  

    
# close the pdf file 
    
if(pflag == 1):           
    pdf.close()              
  