# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:28:00 2015

@author: rogeryelle
"""
from pylab import *
import numpy as np
from matplotlib.pyplot import *
from read_rates import read_rates
from consolidate_rates import consolidate_rates

srun='h'
RT=2575.E5

# read reaction rates

crates = read_rates('../runs/'+srun+'/output/chemrates.out')
prates = read_rates('../runs/'+srun+'/output/photorates.out')
erates = read_rates('../runs/'+srun+'/output/elerates.out')
grates = read_rates('../runs/'+srun+'/output/grates.out')
alt = crates["alt"]
nalt = size(alt)

title = prates["title"] + erates["title"] + grates["title"] + crates["title"] 
crct = np.asarray(crates["rct"], dtype=np.float)
prct = np.asarray(prates["rct"], dtype=np.float)
erct = np.asarray(erates["rct"], dtype=np.float)
grct = np.asarray(grates["rct"], dtype=np.float)
rct = np.concatenate((prct,erct,grct,crct),axis=1)
nrct = size(title)

# calculate column integrated rates

colrates = column_rates(alt,rates)   

# consolidate similar reactions

consolidate_rates(rct,colrates,title)