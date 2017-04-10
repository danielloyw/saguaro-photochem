# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 11:25:28 2014

@author: rogeryelle
"""


import numpy as np

def read_atm(filename):
    
    f=open(filename,'r')
        
    nalt, nmol = map(int,f.readline().split())
#    cos_sza = np.float(f.readline().strip())

    nmol_lines = nmol/10 + 1
    nalt_lines = nalt/10 + 1
    
    line=f.readline()    
    sline='Ntot'
    for x in xrange(0,nmol_lines):
        sline=sline+f.readline()
    name=sline.split()
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        line=f.readline()
        sline=sline+line 
    sline=np.array(sline.split())
    alt=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline()
    sline=np.array(sline.split())
    rad=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    grv=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    Tn=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    Te=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    prs=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    rho=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    mmw=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    Kzz=sline.astype(np.float)
      
    den=np.zeros((nalt,nmol+1), dtype=np.float)
    for n in xrange(0,nmol+1):
        line=f.readline()
        sline=''
        for x in xrange(0,nalt_lines):
            sline=sline+f.readline() 
        sline=np.array(sline.split())
        den[:,n]=sline.astype(np.float)

#    atm = {"nalt":nalt, "nmol":nmol, "cos_sza":cos_sza,"name":name,   \
#        "alt":alt, "rad":rad, "grv":grv, "Tn":Tn, "Te":Te, "prs":prs, \
#        "rho":rho, "mmw":mmw, "Kzz":Kzz, "den":den}
    atm = {"nalt":nalt, "nmol":nmol, "name":name,   \
        "alt":alt, "rad":rad, "grv":grv, "Tn":Tn, "Te":Te, "prs":prs, \
        "rho":rho, "mmw":mmw, "Kzz":Kzz, "den":den}
    atm["nalt"] = nalt
    atm["nmol"] = nmol
#    atm["cos_sza"] = cos_sza
    atm["name"] = name
    atm["alt"] = alt
    atm["rad"] = rad
    atm["Tn"] = Tn
    atm["Te"] = Te
    atm["prs"] = prs
    atm["rho"] = rho
    atm["mmw"] = mmw
    atm["Kzz"] = Kzz
    atm["den"] = den
   
    f.close()
    return atm
