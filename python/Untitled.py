# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:02:09 2015

@author: rogeryelle
"""

f=open('../runs/g/output/atm1D.out','r')
    
nalt, nmol = map(int,f.readline().split())
cos_sza = np.float(f.readline().strip())

nmol_lines = nmol/10 + 1
nalt_lines = nalt/10 + 1

line=f.readline()    
sline='Ntot'
for x in xrange(0,nmol_lines):
    line = f.readline()
    print line
    sline=sline+line       
name=sline.split()

line=f.readline()
sline=''
for x in xrange(0,nalt_lines):
    line=f.readline()
    print line
    sline=sline+line 
sline=np.array(sline.split())
alt=sline.astype(np.float)

line=f.readline()
sline=''
for x in xrange(0,nalt_lines):
    sline=sline+f.readline()
sline=np.array(sline.split())
rad=sline.astype(np.float)
    
