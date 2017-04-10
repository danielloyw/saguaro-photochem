# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 11:25:28 2014

@author: rogeryelle
"""


import numpy as np

def read_ext(filename,smolec):
    
    f=open(filename,'r')
    nprd, nalt = map(int,f.readline().split())
    array=np.zeros((nalt), dtype=np.float)
    ext_prd=np.zeros((nalt), dtype=np.float)
    nalt_lines = nalt/10 + 1
    for n in xrange(0,nprd):
        sline=f.readline()
        sname=sline.strip()
        slist=''
        for x in xrange(0,nalt_lines):
            slist=slist+f.readline()
        slist=np.array(slist.split())
        array=slist.astype(np.float)
        if(sname == smolec):
            ext_prd = array
    
    f.close()
    
    return ext_prd