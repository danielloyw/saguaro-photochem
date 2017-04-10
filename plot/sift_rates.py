# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:02:36 2015

@author: rogeryelle
"""
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from read_rates import read_rates
from read_summary import read_summary

#def sift_rates(srun,smolec):
srun = 'h'
smolec = 'HCN'

# read reaction rates

crates = read_rates('../runs/'+srun+'/output/chemrates.out')
prates = read_rates('../runs/'+srun+'/output/photorates.out')
erates = read_rates('../runs/'+srun+'/output/elerates.out')
grates = read_rates('../runs/'+srun+'/output/grates.out')
title = prates["title"] + erates["title"] + grates["title"] + crates["title"] 
nrct = size(title)
alt = crates["alt"]    

nalt = size(alt)
crct = np.asarray(crates["rct"], dtype=np.float)
prct = np.asarray(prates["rct"], dtype=np.float)
erct = np.asarray(erates["rct"], dtype=np.float)
grct = np.asarray(grates["rct"], dtype=np.float)
rct = np.concatenate((prct,erct,grct,crct),axis=1)
       
# find reactions for smolec

nlss_list = np.zeros(nrct,dtype=int)
nprd_list = np.zeros(nrct,dtype=int)
mp=0
ml=0
for n in xrange(0,nrct):
    nt = title[n].find(smolec)
    if(nt>30):
        nprd_list[mp] = n
        mp += 1 
    elif((nt>-1)&(nt<30)):
        nlss_list[ml] = n
        ml += 1

nprd = mp
nlss = ml        
nprd_list=nprd_list[0:nprd-1]
nlss_list=nlss_list[0:nlss-1]

title_prd = title[nprd_list]
rates_prd = rates[nprd_list]
colrates_prd = colrates[nprd_list]

title_lss = title[nlss_list]
rates_lss = rates[nlss_list]
colrates_lss = colrates[nlss_list]

# consolidate reactions



sr1=[]
sr2=[]
sp1=[]
sp2=[]
sp3=[]
nr=0
for x in title:
    stmp = title[nr].split() 
    sr1.append(stmp[0])
    sr2.append(stmp[2])
    sp1.append(stmp[4])
    sp2.append(stmp[6])
    if (size(stmp)<10):
        sp3.append('')
    else:
        sp3.append(stmp[8])
    nr += 1      

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
#  

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
  k=0
  
  for n in xrange(1, nprd):
    lfind = np.bool_(1)
    m = 0
    while logical_and(np.less_equal(m,k),(lfind)):       
      l1 = (sr1_prd[n] == sr1_new[m])
      l2 = (sr2_prd[n] == sr2_new[m])
      l3 = (sp1_prd[n] == sp1_new[m])
      l4 = (sp2_prd[n] == sp2_new[m])
      l5 = (sp3_prd[n] == sp3_new[m])      
      lchk = np.logical_and(l1,l2)
      lchk = np.logical_and(lchk,l3)
      if logical_and((sp2[n] != ''),(sp3[m] != '-')):
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

  nprd = k + 1
  title_prd = title_new[0:nprd]
  rates_prd = rates_new[0:nalt,0:nprd]
  colrates_prd = colrates_new[0:nprd]
  
  print 'Production Reactions'
  for n in xrange(0,min([5,nprd])):
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
      if logical_and((sp2[n] != ''),(sp3[m] != '-')):
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

