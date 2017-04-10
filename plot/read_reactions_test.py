# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 08:37:15 2014

@author: rogeryelle
"""

import numpy as np
from csv_to_list import csv_to_list

filename = '../runs/g/input/nreactions.csv'

rtab  = csv_to_list(filename)

line0 = rtab[0]
nrct=np.int(line0[0])
rcnt1 = []
rcnt2 = []
prod1 = []
prod2 = []
prod3 = []
itype=np.zeros(nrct, dtype=int)
rk=np.zeros((nrct,10), dtype=float)

for m in range(nrct):
    print m
    n = m + 3
    line = rtab[n]
    rcnt1.append(line[1])
    rcnt2.append(line[2])
    prod1.append(line[3])
    prod2.append(line[4])
    prod3.append(line[5])
    itype[m] = np.int(line[6])
    rk[m,0] = np.float(line[7])
    rk[m,1] = np.float(line[8])
    rk[m,2] = np.float(line[9])
    rk[m,3] = np.float(line[10])
    rk[m,4] = np.float(line[11])
    rk[m,5] = np.float(line[12])
    rk[m,6] = np.float(line[13])
    rk[m,7] = np.float(line[14])
    rk[m,8] = np.float(line[15])
    rk[m,9] = np.float(line[16])

react={"itype":itype,"rcnt1":rcnt1,"rcnt2":rcnt2,"prod1":prod1, \
    "prod2":prod2,"prod3":prod3,"rk":rk}
react["itype"]=itype
react["rcnt1"]=rcnt1
react["rcnt2"]=rcnt2
react["prod1"]=prod1
react["prod2"]=prod2
react["prod3"]=prod3
react["rk"]=rk
