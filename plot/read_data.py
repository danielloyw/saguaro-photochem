# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 18:51:33 2015

@author: rogeryelle
"""
import numpy as np
from operator import sub, add

def read_data(smol):

    f=open('../data/observations/'+smol+'_obs.txt','r')
    
    # number of CIRS observations
    
    sline=f.readline()
    sline.strip()
    sarray=sline.split()
    ncirs=np.int(sarray[1])
    
    # number of UVIS observations
    
    sline=f.readline()
    sline.strip()
    sarray=sline.split()
    nuvis = np.int(sarray[1])
    
    # number of INMS observations
    
    sline=f.readline()
    sline.strip()
    sarray=sline.split()
    ninms=np.int(sarray[1])
    
    # Loop over CIRS observations
    
    ndat_cirs = []
    alt_cirs = []
    tmp_cirs = []
    prs_cirs = []
    xmol_cirs = []
    xmin_cirs = []
    xmax_cirs = []
    for n in xrange(0,ncirs):
    
        for m in xrange(0,8):
            f.readline()
        
        sline=f.readline()
        sline.strip()
        sarray=sline.split()
        ndat = np.int(sarray[0])
        
        f.readline()
        alt = []
        prs = []
        tmp = []
        xmol = []
        xmin = []
        xmax = []
        for m in xrange(0,ndat):
            sline=f.readline()
            sline.strip()
            sarray=sline.split()
            alt.append(np.float(sarray[0]))
            prs.append(np.float(sarray[1]))
            tmp.append(np.float(sarray[2]))
            xmol.append(np.float(sarray[3]))
            xmin.append(np.float(sarray[4]))
            xmax.append(np.float(sarray[5]))
            
        ndat_cirs.append(ndat)
        alt_cirs.append(alt)    
        tmp_cirs.append(tmp)
        prs_cirs.append(prs)
        xmol_cirs.append(xmol)
        xmin_cirs.append(xmin)
        xmax_cirs.append(xmax)
        
    cirs = {"ndat":ndat_cirs, "alt":alt_cirs,"xmol":xmol_cirs,"xmin":xmin_cirs, \
        "xmax":xmax_cirs}        
    cirs["alt"] = alt_cirs
    cirs["xmol"] = xmol_cirs
    cirs["xmin"] = xmin_cirs
    cirs["xmax"] = xmax_cirs
    
    # Loop over UVIS observations
    
    ndat_uvis = []
    alt_uvis = []
    xmol_uvis = []
    xmin_uvis = []
    xmax_uvis = []
    for n in xrange(0,nuvis):
    
        for m in xrange(0,3):
            f.readline()
        
        sline=f.readline()
        sline.strip()
        sarray=sline.split()
        ndat = np.int(sarray[0])
        
        alt = []
        xmol = []
        xerr = []
        for m in xrange(0,ndat):
            sline=f.readline()
            sline.strip()
            sarray=sline.split()
            alt.append(np.float(sarray[0]))
            xmol.append(np.float(sarray[1]))
            xerr.append(np.float(sarray[2]))

        xmin = map(sub,xmol,xerr)
        xmax = map(sub,xmol,xerr)            
        ndat_uvis.append(ndat)
        alt_uvis.append(alt)    
        xmol_uvis.append(xmol)
        xmin_uvis.append(xmin)
        xmax_uvis.append(xmax)

    uvis = {"ndat":ndat_uvis, "alt":alt_uvis,"xmol":xmol_uvis,"xmin":xmin_uvis, \
        "xmax":xmax_uvis}        
    uvis["alt"] = alt_uvis
    uvis["xmol"] = xmol_uvis
    uvis["xmin"] = xmin_uvis
    uvis["xmax"] = xmax_uvis
    
    # Loop over INMS observations
    
    ndat_inms = []
    alt_inms = []
    xmol_inms = []
    xmin_inms = []
    xmax_inms = []
    for n in xrange(0,ninms):
    
        for m in xrange(0,3):
            f.readline()
        
        sline=f.readline()
        sline.strip()
        sarray=sline.split()
        ndat = np.int(sarray[0])
        
        alt = []
        xmol = []
        xerr = []
        for m in xrange(0,ndat):
            sline=f.readline()
            sline.strip()
            sarray=sline.split()
            alt.append(np.float(sarray[0]))
            xmol.append(np.float(sarray[1]))
            xerr.append(np.float(sarray[2]))

        xmin = map(sub,xmol,xerr)
        xmax = map(add,xmol,xerr)            
        ndat_inms.append(ndat)
        alt_inms.append(alt)    
        xmol_inms.append(xmol)
        xmin_inms.append(xmin)
        xmax_inms.append(xmax)

    inms = {"ndat":ndat_inms, "alt":alt_inms,"xmol":xmol_inms,"xmin":xmin_inms, \
        "xmax":xmax_inms}        
    inms["alt"] = alt_inms
    inms["xmol"] = xmol_inms
    inms["xmin"] = xmin_inms
    inms["xmax"] = xmax_inms
        
    f.close()

    data ={"cirs":cirs, "uvis":uvis, "inms":inms}
    
    return data