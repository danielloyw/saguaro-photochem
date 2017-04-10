# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:58:54 2015

@author: rogeryelle
"""
import numpy as np
from plot_molecule import plot_molecule

srun =   raw_input(' Enter run id ..... > ')
srun.strip()

smolec = raw_input(' Enter molecule ... > ')
smolec.strip()

sopt =   raw_input(' 0=screen,1=pdf ... > ')
sopt.strip()
iopt = np.int(sopt)

plot_molecule(srun,smolec,iopt)
