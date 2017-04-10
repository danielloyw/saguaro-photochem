# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:48:51 2015

@author: rogeryelle
"""

from pylab import *
import numpy as np

def column_rates(alt,rct):

    RT=2575.E5
    rz = RT + alt
    rz2 = np.power(rz/RT,2)
    nrct = size(rates,1)
    colrates = np.zeros(nrct, dtype=float)
    nr = 0
    for x in xrange(0,nrct):
       sm = 0.
       nz = 1
       for z in xrange(1,nalt):
           sm = sm + 0.5*(rct[nz,nr]*rz2[nz]+ \
           rct[nz-1,nr]*rz2[nz-1])*(rz[nz]-rz[nz-1])
           nz += 1
       colrates[nr] = sm
       nr += 1
    return colrates
