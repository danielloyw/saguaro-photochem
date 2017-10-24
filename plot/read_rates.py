# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 11:25:28 2014

@author: rogeryelle
"""

import math
import numpy as np

def read_rates(filename):
    title=[]  
    f=open(filename,'r')
        
    nrct, nalt = map(int,f.readline().split())
    rct=np.zeros((nalt,nrct), dtype=np.float)
    nalt_lines = int(math.ceil(float(nalt)/10.))
    f.readline()
    sline=''
    for x in range(0,nalt_lines):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    alt=sline.astype(np.float)
    for n in range(0,nrct):
        sline = f.readline()
        sline.strip()
        title.append(sline.strip())
        slist=''
        for x in range(0,nalt_lines):
            slist=slist+f.readline()
        slist=np.array(slist.split())
        rct[:,n]=slist.astype(np.float)
    
    f.close()
    
    rates = {"nalt":nalt, "nrct":nrct,"title":title,   \
        "alt":alt, "rct":rct}
    rates["nalt"] = nalt
    rates["nrct"] = nrct
    rates["title"] = title
    rates["rct"] = rct
    
    return rates