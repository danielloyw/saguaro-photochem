# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 11:25:28 2014

@author: rogeryelle
"""

import numpy as np    

def read_summary(filename):
        
    alt=np.zeros(301, dtype=float)
    den=np.zeros(301, dtype=float)
    mol=np.zeros(301, dtype=float)
    flx=np.zeros(301, dtype=float)
    pr_ext=np.zeros(301, dtype=float)
    pr_ph=np.zeros(301, dtype=float)
    ls_ph=np.zeros(301, dtype=float)
    pr_pe=np.zeros(301, dtype=float)
    ls_pe=np.zeros(301, dtype=float)
    pr_chem=np.zeros(301, dtype=float)
    ls_chem=np.zeros(301, dtype=float)
    pr_net=np.zeros(301, dtype=float)
    ls_net=np.zeros(301, dtype=float)
    ls_con=np.zeros(301, dtype=float)
    cvg_flx=np.zeros(301, dtype=float)
    bal=np.zeros(301, dtype=float)
    
    f=open(filename,'r')
    
    for x in range(0,14):
        f.readline()
        
    n = 0
    for line in f:
        line = line.split()
        alt[n]=np.float(line[0])
        den[n]=np.float(line[1])
        mol[n]=np.float(line[2])
        flx[n]=np.float(line[3])
        pr_ext[n]=np.float(line[4])
        pr_ph[n]=np.float(line[5])
        ls_ph[n]=np.float(line[6])
        pr_pe[n]=np.float(line[7])
        ls_pe[n]=np.float(line[8])
        pr_chem[n]=np.float(line[9])
        ls_chem[n]=np.float(line[10])
        pr_net[n]=np.float(line[11])
        ls_net[n]=np.float(line[12])
        ls_con[n]=np.float(line[13])
        cvg_flx[n]=np.float(line[14])
        bal[n]=np.float(line[15])       
        n += 1
            
    f.close()
    
    molsum={"alt":alt,"den":den,"mol":mol,"flx":flx,"pr_ext":pr_ext,"pr_ph":pr_ph, \
        "ls_ph":ls_ph,"pr_pe":pr_pe,"ls_pe":ls_pe,"pr_chem":pr_chem, \
        "ls_chem":ls_chem,"pr_net":pr_net,"ls_net":ls_net,"ls_con":ls_con, \
        "cvg_flx":cvg_flx,"bal":bal}
    
    return molsum          