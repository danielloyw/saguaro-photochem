# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:04:34 2015

@author: rogeryelle
"""
from pylab import *
import numpy as np
from matplotlib.pyplot import *

def consolidate_rates(rates,colrates,title):

    nrct = size(title)
    
  #  get reactant and product molecules 

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
        
    lchk = np.bool_(1)
    title_new = title
    sr1_new = sr1
    sr2_new = sr2
    sp1_new = sp1
    sp2_new = sp2
    sp3_new = sp3
    rates_new = rates
    colrates_new = colrates
    k=0
      
    for n in xrange(1, nrct):
        lfind = np.bool_(1)
        m = 0
        while logical_and(np.less_equal(m,k),(lfind)):       
          l1 = (sr1[n] == sr1_new[m])
          l2 = (sr2[n] == sr2_new[m])
          l3 = (sp1[n] == sp1_new[m])
          l4 = (sp2[n] == sp2_new[m])
          l5 = (sp3[n] == sp3_new[m])      
          lchk = np.logical_and(l1,l2)
          lchk = np.logical_and(lchk,l3)
          if logical_and((sp2[n] != ''),(sp3[m] != '-')):
              lchk = logical_and(lchk,l4)
          if logical_and((sp3[n] != ''),(sp3[m] != '-')):
              lchk = logical_and(lchk,l5)        
          if(lchk):
              lfind = np.bool_(0)
              rates_new[:,m] = rates_new[:,m] + rates[:,n]
              colrates_new[m] = colrates_new[m] + colrates[n] 
        
          m += 1   
              
        if(lfind):
            k += 1
            title_new[k] = title[n]
            sr1_new[k] = sr1[n]
            sr2_new[k] = sr2[n]
            sp1_new[k] = sp1[n]
            sp2_new[k] = sp2[n]
            sp3_new[k] = sp3[n]
            rates_new[:,k] = rates[:,n]
            colrates_new[k] = colrates[n]
    
    nrct = k + 1
    title = title_new[0:nrct]
    rates = rates_new[:,0:nrct]
    colrates = colrates_new[0:nrct]
        
    return 