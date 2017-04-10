# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:58:54 2015

@author: rogeryelle
"""
#import numpy as np
from read_molecules import read_molecules
from plot_molecule import plot_molecule

srun = raw_input(' Enter run id > ')
srun.strip()

# read list of molecules

nmolec = read_molecules('../runs/'+srun+'/input/nmolecules.dat') 
imolec = read_molecules('../runs/'+srun+'/input/imolecules.dat') 
name = nmolec["name"]
name.extend(imolec["name"])
stat = nmolec["stat"]
stat.extend(imolec["stat"])

n = 0
for x in name:
    if(stat[n] > 0):
        print name[n]
        plot_molecule(srun,name[n],1)
    n += 1
