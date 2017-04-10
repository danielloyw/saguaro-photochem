# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:28:00 2015

@author: rogeryelle
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from read_rates import read_rates
from read_summary import read_summary

def plot_molecule_new(srun,smolec,pflag):

    RT = 2575.E5
    colrates_lim = 1.E3

    # read molecule summary files
    
    if(os.path.isfile('../runs/'+srun+'/output/molecules/'+smolec+'.OUT')):
        molsum = read_summary('../runs/'+srun+'/output/molecules/'+smolec+'.OUT')   
        ierr = 0
    else:
        ierr = 1
        return

    #  Define PDF filename     

    if(pflag == 1):
        sfileplot = '../runs/'+srun+'/plots/'+smolec+'.pdf'
    
    # read reaction rates
    
    crates = read_rates('../runs/'+srun+'/output/chemrates.out')
    prates = read_rates('../runs/'+srun+'/output/photorates.out')
    erates = read_rates('../runs/'+srun+'/output/elerates.out')
    grates = read_rates('../runs/'+srun+'/output/grates.out')
    alt = crates["alt"]
    nalt = size(alt)
    
    title = prates["title"] + erates["title"] + grates["title"] + crates["title"] 
    nrct = size(title)
    
    crct = np.asarray(crates["rct"], dtype=np.float)
    prct = np.asarray(prates["rct"], dtype=np.float)
    erct = np.asarray(erates["rct"], dtype=np.float)
    grct = np.asarray(grates["rct"], dtype=np.float)
    
    rct = np.concatenate((prct,erct,grct,crct),axis=1)
    
    # find reactions involving smolec
    
    indx = sift_rates(title,smolec)
    
    # create separate arrays for production and loss
    
    title_prd = title(indx[:,0])
    rct_prd = rct[:,indx[:,0]]
    title_lss = title(indx[:,1])
    rct_lss = rct[:,indx[:,1]]
    
    # calculate column integrated rates
    
    colrates_prd = column_rates(alt,rct_prd)
    colrates_lss = column_rates(alt,rct_lss)
    
    # consolidate similar reactions
    
    consolidate_rates(rct_prd,colrates_prd,title_prd)
    consolidate_rates(rct_lss,colrates_lss,title_lss)
    
    # sort reactions in order of decreasing column rate
    
    sort_rates(rct_prd,colrates_prd,title_prd)
    sort_rates(rct_lss,colrates_lss,title_lss)
       
      
    print 'Production Reactions'
    for n in xrange(0,min([5,nprd])):
        print("{0:10.2e}:  {1:s}".format(colrates_prd[n],title_prd[n]))
    
    print 'Loss Reactions'
    for n in xrange(0,min([5,nlss])):
        print("{0:10.2e}:  {1:s}".format(colrates_lss[n],title_lss[n]))
    
    #             Make Plots
    
    alt = 1.E-5*alt
    
    plt.rcParams['xtick.labelsize'] = 8 
    plt.rcParams['ytick.labelsize'] = 8     
    clrs = ['red','blue','green','magenta','cyan']
    ymin = 0
    ymax = 1500
    
    #  Open PDF file
 
    if(pflag == 1):       
        pdf = PdfPages(sfileplot)
    
    fig = plt.figure()    
    
    # plot density upper left
    
    ax1 = fig.add_subplot(2, 2, 1)
    xmax = math.pow(10.,int(math.log10(np.amax(molsum["den"])))+1)
    xmin = xmax/1.E10
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    ax1.semilogx(molsum["den"],molsum["alt"])
    ax1.set_ylabel(r'Altitude (Km)',fontsize='8')
    ax1.set_xlabel(r'Density (cm$^{-3}$)',fontsize='8')          
    ax1.text(xmax,ymax+500,smolec,fontsize='10')                
        
    # plot mole fraction upper right
    
    ax2 = fig.add_subplot(2, 2, 2)
    xmax = math.pow(10.,int(math.log10(np.amax(molsum["mol"])))+1)
    xmin = xmax/1.E10
    ax2.set_xlim(xmin,xmax)
    ax2.set_ylim(ymin,ymax)
    ax2.semilogx(molsum["mol"],molsum["alt"])
    ax2.set_ylabel(r'Altitude (Km)',fontsize='8')
    ax2.set_xlabel(r'Mole Fraction (V/V)',fontsize='8')                          
    
    # plot prodction and loss, Lower left
    
    pmax = np.amax(molsum["pr_net"])
    lmax = np.amax(molsum["ls_net"])
    cmax = np.amax(abs(molsum["cvg_flx"]))
    xmax = np.amax([pmax,lmax,cmax])
    if(xmax > 1.E-10):
        xmax = math.pow(10.,int(math.log10(xmax))+1)
        xmin = xmax/1.E10        
        ax3 = fig.add_subplot(2, 2, 3)
        ax3.set_xlim(xmin,xmax)
        ax3.set_ylim(ymin,ymax)
        ax3.semilogx(molsum["pr_net"],molsum["alt"],color=clrs[0])
        ax3.semilogx(molsum["ls_net"],molsum["alt"],color=clrs[1])
        ax3.semilogx(molsum["cvg_flx"],molsum["alt"],color=clrs[2])
        ax3.semilogx(-1*molsum["cvg_flx"],molsum["alt"],color=clrs[2],linestyle=':')
        ax3.set_ylabel(r'Altitude (Km)',fontsize='8')
        ax3.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$',fontsize='8')          
        ax3.semilogx([3*xmin,30*xmin],[ymax-100,ymax-100],color=clrs[0])
        ax3.text(50*xmin,ymax-100,r'Net Chemical Production',color=clrs[0],fontsize='6')                
        ax3.semilogx([3*xmin,30*xmin],[ymax-160,ymax-160],color=clrs[1])
        ax3.text(50*xmin,ymax-160,r'Net Chemical Loss',color=clrs[1],fontsize='6')       
        ax3.semilogx([3*xmin,30*xmin],[ymax-220,ymax-220],color=clrs[2])         
        ax3.text(50*xmin,ymax-220,r'Flux Divergence',color=clrs[2],fontsize='6')                
        ax3.semilogx([3*xmin,30*xmin],[ymax-280,ymax-280],color=clrs[2],linestyle=':')         
        ax3.text(50*xmin,ymax-280,r'Flux Convergence',color=clrs[2],fontsize='6')                
    
    # plot Fluxes, Lower right
    
    altflx = molsum["alt"]
    flx = molsum["flx"]
    if(np.amax(abs(flx))>1.E-10):
        upalt = altflx[(flx>0)]
        upflx = flx[(flx>0)]
        dnalt = altflx[(flx<0)]
        dnflx = -1*flx[(flx<0)]
        umax = np.amax(upflx)
        dmax = np.amax(dnflx)
        xmax = np.amax([umax,dmax])
        xmax = math.pow(10,int(math.log10(xmax))+1)
        xmin = xmax/1.E10
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.set_xlim(xmin,xmax)
        ax4.set_ylim(ymin,ymax)
        ax4.semilogx(upflx,upalt,color=clrs[0])
        ax4.semilogx(dnflx,dnalt,color=clrs[0],linestyle=':')
        ax4.set_ylabel(r'Altitude (Km)',fontsize='8')
        ax4.set_xlabel(r'Flux (cm$^{-2}$s$^{-1}$',fontsize='8')   
        ax4.semilogx([3*xmin,30*xmin],[ymax-100,ymax-100],color=clrs[0])         
        ax4.text(50*xmin,ymax-100,r'Upward Flux',color=clrs[0],fontsize='6')                
        ax4.semilogx([3*xmin,30*xmin],[ymax-160,ymax-160],color=clrs[0],linestyle=':')         
        ax4.text(50*xmin,ymax-160,r'Downward Flux',color=clrs[0],fontsize='6')                
                       

    if(pflag == 1):
        pdf.savefig(fig)
        plt.close()
    else:
        plt.show()  
   
    # plot production reactions
    
    l = 0
    if(nprd>0):

        fig = plt.figure()      
        
        # upper left

        l1 = l 
        l2 = max([l1,min([l + 5,nprd-1])])
        if (l1 == l2):
            xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1])))+1)
        else:    
            xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1:l2+1])))+1)
        xmin = xmax/1.E10
        
        ax1 = fig.add_subplot(2, 2, 1)
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)
        ax1.set_ylabel(r'Altitude (Km)',fontsize='8')
       
        k = 0
        while ((k<5) & (l<nprd)):
            ax1.semilogx(rates_prd[:,l],alt,color=clrs[k])
            ax1.text(3*xmin,ymax-100-k*60,title_prd[l],color=clrs[k],fontsize='6')
            l = l + 1
            k = k + 1
            
        if (l == nprd): 
            ax1.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')                          
    
        # upper right
        
        if(l<nprd):
    
            l1 = l 
            l2 = max([l1,min([l + 5,nprd-1])])
            if (l1 == l2):
                xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1])))+1)
            else:    
                xmax = np.power(10.,int(np.log10(np.amax(rates_prd[:,l1:l2+1])))+1)
            xmin = xmax/1.E10
            ax2 = fig.add_subplot(2, 2, 2)
            ax2.set_xlim(xmin,xmax)
            ax2.set_ylim(ymin,ymax)
            k = 0
            while ((k<5) & (l<nprd)):
                ax2.semilogx(rates_prd[:,l],alt,color=clrs[k])
                ax2.text(3*xmin,ymax-100-k*60,title_prd[l],color=clrs[k],fontsize='6')
                l = l + 1
                k = k + 1
                
            if (l == nprd): 
                ax1.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')                 
                ax2.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')                      
        
        # lower left
        
        if(l<nprd):
    
            l1 = l 
            l2 = max([l1,min([l + 5,nprd-1])])
            if (l1 == l2):
                xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1])))+1)
            else:    
                xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1:l2+1])))+1)
            xmin = xmax/1.E10
            ax3 = fig.add_subplot(2, 2, 3)
            ax3.set_xlim(xmin,xmax)
            ax3.set_ylim(ymin,ymax)
            ax3.set_ylabel(r'Altitude (Km)',fontsize='8')
            ax3.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')
            k = 0
            while ((k<5) & (l<nprd)):
                ax3.semilogx(rates_prd[:,l],alt,color=clrs[k])
                ax3.text(3*xmin,ymax-100-k*60,title_prd[l],color=clrs[k],fontsize='6')
                l = l + 1
                k = k + 1
                
            if (l == nprd):
                ax2.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')            
        
        # lower right
        
        if(l<nprd):
     
            l1 = l
            l2 = max([l1,min([l + 5,nprd-1])])
            if (l1 == l2):
                xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1])))+1)
            else:    
                xmax = math.pow(10.,int(math.log10(np.amax(rates_prd[:,l1:l2+1])))+1)
            xmin = xmax/1.E10
            ax4 = fig.add_subplot(2, 2, 4)
            ax4.set_xlim(xmin,xmax)
            ax4.set_ylim(ymin,ymax)
            ax4.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')
            k = 0
            while ((k<5) & (l<nprd)):
                ax4.semilogx(rates_prd[:,l],alt,color=clrs[k])
                ax4.text(3*xmin,ymax-100-k*60,title_prd[l],color=clrs[k],fontsize='6')
                l = l + 1
                k = k + 1

        if(pflag == 1):
            pdf.savefig(fig)
            plt.close()
        else:
            plt.show()  

    #
    #   Plot Loss Rates
    #
        
    
    if(nlss>0):
                
        l = 0
        fig = plt.figure()  # create a figure object
        
        # Upper left
        
        l1 = l 
        l2 = max([l1,min([l + 5,nlss-1])])
        if (l1 == l2):
            xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1])))+1)
        else:    
            xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1:l2+1])))+1)
        xmin = xmax/1.E10
        
        ax1 = fig.add_subplot(2, 2, 1)
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)

        k = 0
        while ((k<5) & (l<nlss)):
            ax1.semilogx(rates_lss[:,l],alt,color=clrs[k])
            ax1.text(3*xmin,ymax-100-k*60,title_lss[l],color=clrs[k],fontsize='6')
            l = l + 1
            k = k + 1
            
        ax1.set_ylabel(r'Altitude (Km)',fontsize='8')
        if (l == nlss): 
            ax1.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')   
              
        # Upper Right

        if(l<nlss):
            l1 = l 
            l2 = max([l1,min([l + 5,nlss-1])])
            if (l1 == l2):
                xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1])))+1)
            else:                 
                xmax = np.power(10.,int(np.log10(np.amax(rates_lss[:,l1:l2+1])))+1)
            xmax = max([xmax,1.E-20])
            xmin = xmax/1.E10
            
            ax2 = fig.add_subplot(2, 2, 2)
            ax2.set_xlim(xmin,xmax)
            ax2.set_ylim(ymin,ymax)
            k = 0
            while ((k<5) & (l<nlss)):
                ax2.semilogx(rates_lss[:,l],alt,color=clrs[k])
                ax2.text(3*xmin,ymax-100-k*60,title_lss[l],color=clrs[k],fontsize='6')
                l = l + 1
                k = k + 1
                
            if (l == nlss): 
                ax1.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')                 
                ax2.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')            
        
        
        # lower left
    
        if(l<nlss):    
            l1 = l 
            l2 = max([l1,min([l + 5,nlss-1])])
            if (l1 == l2):
                xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1])))+1)
            else:    
                xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1:l2+1])))+1)
            xmin = xmax/1.E10

            ax3 = fig.add_subplot(2, 2, 3)
            ax3.set_xlim(xmin,xmax)
            ax3.set_ylim(ymin,ymax)
            k = 0
            while ((k<5) & (l<nlss)):
                ax3.semilogx(rates_lss[:,l],alt,color=clrs[k])
                ax3.text(3*xmin,ymax-100-k*60,title_lss[l],color=clrs[k],fontsize='6')
                l = l + 1
                k = k + 1
            
            ax3.set_ylabel(r'Altitude (Km)',fontsize='8')
            ax3.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')
            if (l == nlss):
                ax2.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8') 
        
        # lower right
        
        if(l<nlss):
            l1 = l 
            l2 = max([l1,min([l + 5,nlss-1])])
            if (l1 == l2):
                xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1])))+1)
            else:    
                xmax = math.pow(10.,int(math.log10(np.amax(rates_lss[:,l1:l2+1])))+1)
            xmin = xmax/1.E10
            
            ax4 = fig.add_subplot(2, 2, 4)
            ax4.set_xlim(xmin,xmax)
            ax4.set_ylim(ymin,ymax)
            k = 0
            while ((k<5) & (l<nlss)):
                ax4.semilogx(rates_lss[:,l],alt,color=clrs[k])
                ax4.text(3*xmin,ymax-100-k*60,title_lss[l],color=clrs[k],fontsize='6')
                l = l + 1
                k = k + 1
            
            ax4.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)',fontsize='8')

        if(pflag == 1):
            pdf.savefig(fig)
            plt.close()
        else:
            plt.show()  
    
    # close the pdf file for this molecule
    if(pflag == 1):           
        pdf.close()              
        
    return ierr