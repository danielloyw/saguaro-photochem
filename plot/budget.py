# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:58:54 2015

@author: rogeryelle
"""
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from read_molecules import read_molecules
from read_rates import read_rates
from read_summary import read_summary

R = 3390
colrates_lim = 1.E-3

plt.rcParams['xtick.labelsize'] = 8 
plt.rcParams['ytick.labelsize'] = 8 
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['figure.titlesize'] = 10
plt.rcParams['font.size'] = 4

srun = raw_input(' Enter run id > ')
srun.strip()

smolec = raw_input(' Enter molecule > ')
smolec.strip()

sfileplot = '../runs/'+srun+'/plots/'+smolec+'.pdf'

# read list of molecules

nmolec = read_molecules('../runs/'+srun+'/input/nmolecules.dat') 
imolec = read_molecules('../runs/'+srun+'/input/imolecules.dat') 
molec = [nmolec,imolec]

# read reaction rates

crates = read_rates('../runs/'+srun+'/output/chemrates.out')
prates = read_rates('../runs/'+srun+'/output/photorates.out')
erates = read_rates('../runs/'+srun+'/output/elerates.out')
alt = crates["alt"]
nalt = size(alt)

title = prates["title"] + erates["title"] + crates["title"] 
nrct = size(title)

crct = np.asarray(crates["rct"], dtype=np.float)
prct = np.asarray(prates["rct"], dtype=np.float)
erct = np.asarray(erates["rct"], dtype=np.float)

rct = np.concatenate((prct,erct,crct),axis=1)

# calculate column integrated rates

alt = alt/1.E5
rz = R + alt
rz2 = np.power(rz/R,2)
colrates = np.zeros(nrct, dtype=float)
for nr in xrange(0,nrct):
   sm = 0.
   for nz in xrange(1,nalt):
       sm = sm + 0.5*(rct[nz,nr]*rz2[nz]+ \
       rct[nz-1,nr]*rz2[nz-1])*(rz[nz]-rz[nz-1])
   colrates[nr] = sm

# read molecule summary files

molsum = read_summary('../runs/'+srun+'/output/molecules/'+smolec+'.OUT')   

# get reactant and product molecules from reaction title

sr1=[]
sr2=[]
sp1=[]
sp2=[]
sp3=[]
for nr in range(0,nrct):
    stmp = title[nr].split() 
    sr1.append(stmp[0])
    sr2.append(stmp[2])
    sp1.append(stmp[4])
    sp2.append(stmp[6])
    if (size(stmp)<10):
        sp3.append('')
    else:
        sp3.append(stmp[8])

nlss_list = []
nprd_list = []
for nr in range(0,nrct):

   if(colrates[nr]>colrates_lim):
       
       # find reactions with molecule as a reactant

       if((smolec==sr1[nr]) or (smolec==sr2[nr])):
           nlss_list.append(nr)

       # find reactions with molecule as a product

       if((smolec==sp1[nr]) or (smolec==sp2[nr]) or (smolec==sp3[nr])):
           nprd_list.append(nr)


nprd = size(nprd_list)
nlss = size(nlss_list)

title = np.array(title)
sr1 = np.array(sr1)
sr2 = np.array(sr2)
sp1 = np.array(sp1)
sp2 = np.array(sp2)
sp3 = np.array(sp3)
nprd_indx = np.asarray(nprd_list, dtype=int)
nprd_indx.shape = (nprd)
nlss_indx = np.asarray(nlss_list, dtype=int)
nlss_indx.shape = (nlss)

#
#  switch to first consolidate then to sort
#  make selection, sorting, and consolidation into functions

# sort production reactions by rate

if (nprd > 0):
    
  # sort reactions
    
  ftmp_prd = colrates[nprd_indx]
  rtmp_prd = rct[:,nprd_indx]
  stmp_prd = title[nprd_indx]
  stmp_r1 = sr1[nprd_indx]
  stmp_r2 = sr2[nprd_indx]
  stmp_p1 = sp1[nprd_indx]
  stmp_p2 = sp2[nprd_indx]
  stmp_p3 = sp3[nprd_indx]
  sort_indx_prd = flipud(np.argsort(ftmp_prd))
  colrates_prd = ftmp_prd[sort_indx_prd]
  rates_prd = rtmp_prd[:,sort_indx_prd]
  title_prd = stmp_prd[sort_indx_prd]
  sr1_prd = stmp_r1[sort_indx_prd]
  sr2_prd = stmp_r2[sort_indx_prd]
  sp1_prd = stmp_p1[sort_indx_prd]
  sp2_prd = stmp_p2[sort_indx_prd]
  sp3_prd = stmp_p3[sort_indx_prd]
    
  #  consolidate similar reactions 
    
  lchk = np.bool_(1)
  title_new = title_prd
  sr1_new = sr1_prd
  sr2_new = sr2_prd
  sp1_new = sp1_prd
  sp2_new = sp2_prd
  sp3_new = sp3_prd  
  rates_new = rates_prd
  colrates_new = colrates_prd
  k=0 # consolidated array index
  
  for n in range(1, nprd):
    lfind = np.bool_(1) # new reaction?
    m = 0 # reference reaction in consolidated array
    while logical_and(np.less_equal(m,k),(lfind)):
      l1 = (sr1_prd[n] == sr1_new[m])
      l2 = (sr2_prd[n] == sr2_new[m])
      l3 = (sp1_prd[n] == sp1_new[m])
      l4 = (sp2_prd[n] == sp2_new[m])
      l5 = (sp3_prd[n] == sp3_new[m])
      lchk = np.logical_and(l1,l2)
      lchk = np.logical_and(lchk,l3)
      if logical_and((sp2[n] != ''),(sp2[m] != '-')):
          lchk = logical_and(lchk,l4)
      if logical_and((sp3[n] != ''),(sp3[m] != '-')):
          lchk = logical_and(lchk,l5)
      if(lchk):
          lfind = np.bool_(0)
          rates_new[:,m] = rates_new[:,m] + rates_prd[:,n]
          colrates_new[m] = colrates_new[m] + colrates_prd[n] 

      m += 1   
          
    if(lfind):
        k += 1
        title_new[k] = title_prd[n]
        sr1_new[k] = sr1_prd[n]
        sr2_new[k] = sr2_prd[n]
        sp1_new[k] = sp1_prd[n]
        sp2_new[k] = sp2_prd[n]
        sp3_new[k] = sp3_prd[n]
        rates_new[:,k] = rates_prd[:,n]
        colrates_new[k] = colrates_prd[n]

  nprd = k + 1 # size is k+1 because python starts at 0
  title_prd = title_new[0:nprd]
  rates_prd = rates_new[0:nalt,0:nprd]
  colrates_prd = colrates_new[0:nprd]
  
  print 'Production Reactions'
  for n in range(0,min([5,nprd])):
     print("{0:10.2e}:  {1:s}".format(colrates_prd[n],title_prd[n]))

if (nlss > 0):
    
  # sort reactions
    
  ftmp_lss = colrates[nlss_indx]
  rtmp_lss = rct[:,nlss_indx]
  stmp_lss = title[nlss_indx]
  stmp_r1 = sr1[nlss_indx]
  stmp_r2 = sr2[nlss_indx]
  stmp_p1 = sp1[nlss_indx]
  stmp_p2 = sp2[nlss_indx]
  stmp_p3 = sp3[nlss_indx]
  sort_indx_lss = flipud(np.argsort(ftmp_lss))
  colrates_lss = ftmp_lss[sort_indx_lss]
  rates_lss = rtmp_lss[:,sort_indx_lss]
  title_lss = stmp_lss[sort_indx_lss]
  sr1_lss = stmp_r1[sort_indx_lss]
  sr2_lss = stmp_r2[sort_indx_lss]
  sp1_lss = stmp_p1[sort_indx_lss]
  sp2_lss = stmp_p2[sort_indx_lss]
  sp3_lss = stmp_p3[sort_indx_lss]
    
  #  consolidate similar reactions
    
  lchk = np.bool_(1)
  title_new = title_lss
  sr1_new = sr1_lss
  sr2_new = sr2_lss
  sp1_new = sp1_lss
  sp2_new = sp2_lss
  sp3_new = sp3_lss  
  rates_new = rates_lss
  colrates_new = colrates_lss
  k=0
  
  for n in xrange(1, nlss):
      
    lfind = np.bool_(1)
    m = 0
    while logical_and(np.less_equal(m,k),(lfind)):
      l1 = (sr1_lss[n] == sr1_new[m])
      l2 = (sr2_lss[n] == sr2_new[m])
      l3 = (sp1_lss[n] == sp1_new[m])
      l4 = (sp2_lss[n] == sp2_new[m])
      l5 = (sp3_lss[n] == sp3_new[m])      
      lchk = np.logical_and(l1,l2)
      lchk = np.logical_and(lchk,l3)
      if logical_and((sp2[n] != ''),(sp2[m] != '-')):
          lchk = logical_and(lchk,l4)
      if logical_and((sp3[n] != ''),(sp3[m] != '-')):
          lchk = logical_and(lchk,l5)        
      if(lchk):
          lfind = np.bool_(0)
          rates_new[:,m] = rates_new[:,m] + rates_lss[:,n]
          colrates_new[m] = colrates_new[m] + colrates_lss[n] 
      m += 1
          
    if(lfind):
        k += 1
        title_new[k] = title_lss[n]
        sr1_new[k] = sr1_lss[n]
        sr2_new[k] = sr2_lss[n]
        sp1_new[k] = sp1_lss[n]
        sp2_new[k] = sp2_lss[n]
        sp3_new[k] = sp3_lss[n]
        rates_new[:,k] = rates_lss[:,n]
        colrates_new[k] = colrates_lss[n]

  nlss = k + 1
  title_lss = title_new[0:nlss]
  rates_lss = rates_new[0:nalt,0:nlss]
  colrates_lss = colrates_new[0:nlss]
  print 'Loss Reactions'
  for n in xrange(0,min([5,nlss])):
     print("{0:10.2e}:  {1:s}".format(colrates_lss[n],title_lss[n]))

fo = open('../runs/'+srun+'/plots/test.txt', "wb")
fo.write(str(rates_lss[:,0]));
fo.close()

##
##               Make Plots
##

clrs = ['red','blue','green','magenta','cyan']
ymin = 0
ymax = 250

#  Open PDF file

pdf = PdfPages(sfileplot)

fig = plt.figure()

# plot density upper left

ax1 = fig.add_subplot(2, 2, 1)
xmax = math.pow(10,int(math.log10(np.amax(molsum["den"])))+1)
xmin = xmax/1.E10
ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
ax1.semilogx(molsum["den"],molsum["alt"])
ax1.set_ylabel(r'Altitude (km)')
ax1.set_xlabel(r'Density (cm$^{-3}$)')
ax1.xaxis.set_label_position('top')

# plot mole fraction upper right

ax2 = fig.add_subplot(2, 2, 2)
xmax = math.pow(10,int(math.log10(np.amax(molsum["mol"])))+1)
xmin = xmax/1.E8
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin,ymax)
ax2.semilogx(molsum["mol"],molsum["alt"])
ax2.set_xlabel(r'Mole Fraction (V/V)')
ax2.xaxis.set_label_position('top')

# plot production and loss, Lower left

ax3 = fig.add_subplot(2, 2, 3)
pmax = np.amax(molsum["pr_net"])
lmax = np.amax(molsum["ls_net"])
cmax = np.amax(molsum["cvg_flx"])
dmax = np.amax(-1*molsum["cvg_flx"])
xmax = np.amax([pmax,lmax,cmax,dmax])
xmax = math.pow(10,int(math.log10(xmax))+1)
xmin = 1.
ax3.set_xlim(xmin,xmax)
ax3.set_ylim(ymin,ymax)
ax3.semilogx(molsum["pr_net"],molsum["alt"],color=clrs[0])
ax3.semilogx(molsum["ls_net"],molsum["alt"],color=clrs[1])
ax3.semilogx(molsum["cvg_flx"],molsum["alt"],color=clrs[2])
ax3.semilogx(-1*molsum["cvg_flx"],molsum["alt"],color=clrs[3])
ax3.set_ylabel(r'Altitude (Km)')
ax3.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)')
ax3.text(0.8,0.95-0.08,"Net Prod",color=clrs[0],transform=ax3.transAxes)
ax3.text(0.8,0.95-0.16,"Net Loss",color=clrs[1],transform=ax3.transAxes)
ax3.text(0.8,0.95-0.24,"del flx",color=clrs[2],transform=ax3.transAxes)
ax3.text(0.8,0.95-0.32,"-del flx",color=clrs[3],transform=ax3.transAxes)

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
    ax4.semilogx(dnflx,dnalt,color=clrs[1])
    ax4.set_xlabel(r'Flux (cm$^{-2}$s$^{-1}$)')

pdf.savefig(fig)
plt.close()

# plot production reactions

l=0
xmin = 1.E-6
while (l < nprd-1):
    fig = plt.figure()  # create a figure object
    fig.suptitle('Production Rates')
    for subplot_i in range(1,5):
        ax = fig.add_subplot(2, 2, subplot_i)
        l1 = l
        l2 = max([l1,min([l + 5,nprd-1])])
        xmax = math.pow(10,int(math.log10(np.amax(rates_prd[:,l1:l2+1])))+1)
        xmin = min(xmin,xmax/1.E10)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        k = 0
        while ((k<5) and (l<nprd-1)):
            ax.semilogx(rates_prd[:,l],alt,color=clrs[k])
            ax.text(0.1,0.02+k*0.05,title_prd[l],color=clrs[k],transform=ax.transAxes)
            l = l + 1
            k = k + 1
        if (subplot_i == 1) or (subplot_i == 3):
            ax.set_ylabel(r'Altitude (Km)')
        if (subplot_i == 3) or (subplot_i == 4):
            ax.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)')
        if (l == nprd-1): 
            break
    pdf.savefig()
    plt.close()

    # plot loss reactions

l=0
while (l < nlss-1):
    fig = plt.figure()  # create a figure object
    fig.suptitle('Loss Rates')
    for subplot_i in range(1,5):
        ax = fig.add_subplot(2, 2, subplot_i)
        l1 = l
        l2 = max([l1,min([l + 5,nprd-1])])
        xmax = math.pow(10,int(math.log10(np.amax(rates_lss[:,l1:l2+1])))+1)
        xmin = min(xmin,xmax/1.E10)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        k = 0
        while ((k<5) and (l<nlss-1)):
            l = l + 1
            ax.semilogx(rates_lss[:,l],alt,color=clrs[k])
            ax.text(0.1,0.02+k*0.05,title_lss[l],color=clrs[k],transform=ax.transAxes)
            k = k + 1
        if (subplot_i == 1) or (subplot_i == 3):
            ax.set_ylabel(r'Altitude (Km)')
        if (subplot_i == 3) or (subplot_i == 4):
            ax.set_xlabel(r'Rate (cm$^{-3}$s$^{-1}$)')
        if (l == nlss-1): 
            break
    pdf.savefig()
    plt.close()

pdf.close()