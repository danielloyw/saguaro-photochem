#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read various files in /input/ and /output/
"""

__author__ = "Daniel Lo, Roger Yelle"

import numpy as np
import pandas as pd


# Read atm1D.in and atm1D.out
def atm(filename):
    with open(filename,'r') as f:
        
        # read in number of altitude bins and species
        n_alt, n_mol = map(int, f.readline().split())
        
        # number of lines in each altitude/species block
        n_alt_lines = int(np.ceil(n_alt/10.))
        n_mol_lines = int(np.ceil(n_mol/10.))
        
        f.readline() # 'MOLECULES'
        t_block = 'N_total'
        for i_line in range(0, n_mol_lines):
            t_block = t_block + f.readline()
        name = t_block.split()
        
        f.readline() # 'ALTITUDE (cm)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        alt = t_block.astype(float)
        
        f.readline() # 'RADIUS (cm)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        rad = t_block.astype(float)
        
        f.readline() # 'GRAVITY (cm s-2)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        grv = t_block.astype(float)
        
        f.readline() # 'NEUTRAL TEMPERATURE (Kelvins)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        Tn = t_block.astype(float)
        
        f.readline() # 'ELECTRON TEMPERATURE (Kelvins)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        Te = t_block.astype(float)
        
        f.readline() # 'PRESSURE (dyne/cm^2)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        prs = t_block.astype(float)
        
        f.readline() # 'MASS DENSITY (g cm^-3)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        rho = t_block.astype(float)
        
        f.readline() # 'MEAN MOLECULAR WEIGHT (amu)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        mmw = t_block.astype(float)
        
        f.readline() # 'EDDY COEFFICIENT (cm^2s^-1)'
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        Kzz = t_block.astype(float)
          
        den = np.zeros((n_alt, n_mol+1))
        for i_mol in range(0, n_mol+1):
            f.readline() # Species name
            t_block = ''
            for i_line in range(0, n_alt_lines):
                t_block = t_block + f.readline() 
            t_block = np.array(t_block.split())
            den[:, i_mol] = t_block.astype(float)
    
        profiles = pd.DataFrame({'Altitude':alt, 'Radius':rad, 'Gravity':grv, 
                             'Tn':Tn, 'Te':Te, 'P':prs, 'Mass Density': rho, 
                             'MMW':mmw, 'Kzz':Kzz})
        
        for i_mol in range(0, n_mol+1):
            profiles[name[i_mol]] = den[:, i_mol]        
        
        atm = {'n_alt':n_alt, 'n_mol':n_mol, 
               'Species':name, 'Profiles':profiles}
   
    return atm


# Read files for specific species (*.out files in /molecules/)
def summary(filename):
        
    # initialize arrays
    alt = []
    den = []
    mol = []
    flx = []
    pr_ext = []
    pr_ph = []
    ls_ph = []
    pr_pe = []
    ls_pe = []
    pr_chem = []
    ls_chem = []
    pr_net = []
    ls_net = []
    ls_con = []
    cvg_flx = []
    bal = []
    
    with open(filename, 'r') as f:
        
        # skip lines with column values and header labels
        for i in range(0,14):
            f.readline()
        
        # read in altitude profiles
        for line in f:
            line = line.split()
            alt.append(float(line[0]))
            den.append(float(line[1]))
            mol.append(float(line[2]))
            flx.append(float(line[3]))
            pr_ext.append(float(line[4]))
            pr_ph.append(float(line[5]))
            ls_ph.append(float(line[6]))
            pr_pe.append(float(line[7]))
            ls_pe.append(float(line[8]))
            pr_chem.append(float(line[9]))
            ls_chem.append(float(line[10]))
            pr_net.append(float(line[11]))
            ls_net.append(float(line[12]))
            ls_con.append(float(line[13]))
            cvg_flx.append(float(line[14]))
            bal.append(float(line[15]))
    
    return pd.DataFrame({'Altitude':alt, 'Density':den, 'Mixing Ratio':mol, 
                         'Flux':flx, 'Production-External':pr_ext, 
                         'Production (Photo)':pr_ph, 'Loss (Photo)':ls_ph, 
                         'Production (Electron)':pr_pe, 
                         'Loss (Electron)':ls_pe, 
                         'Production (Chemistry)':pr_chem, 
                         'Loss (Chemistry)':ls_chem, 
                         'Production (Net)':pr_net, 'Loss (Net)':ls_net, 
                         'Loss (Condensation)':ls_con, 'cvg_flx':cvg_flx, 
                         'Balance':bal})


# Read nmolecules.dat and imolecules.dat
def molecules(filename):
    with open(filename, 'r') as f:
    
        f.readline() # header for number of species
        f.readline() # column labels
        
        # initialize arrays
        name = []
        stat = []
        chrg = []
        wght = []
        hyd = []
        car = []
        n14 = []
        n15 = []
        oxy = []
    
        for line in f:
            line=line.split()
            name.append(line[1])
            stat.append(int(line[2]))
            chrg.append(int(line[3]))
            wght.append(float(line[4]))
            hyd.append(int(line[5]))
            car.append(int(line[6]))
            n14.append(int(line[7]))
            n15.append(int(line[8]))
            oxy.append(int(line[9]))
        f.close()
        
    return pd.DataFrame({'Species':name, 'stat':stat, 'Charge':chrg, 
                     'Weight':wght, 
                     'nH':hyd, 'nC':car, 'nN14':n14, 'nN15':n15, 'nO':oxy})


# Read chemrates.out, photorates.out, elerates.out, ratecoeff.out
def rates(filename):
    
    # initialize array for reaction names
    title = []  
    
    with open(filename, 'r') as f:
        # read in number of reactions and altitude bins
        n_rct, n_alt = map(int, f.readline().split())
        
        # initialize rate array
        rct = np.zeros((n_alt, n_rct))
        
        # number of lines in each altitude/reaction block
        n_alt_lines = int(np.ceil(n_alt/10.))
        
        f.readline() # altitude header
        
        # read in altitude block
        t_block = ''
        for i_line in range(0, n_alt_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        alt = t_block.astype(float)
        
        # read in reaction blocks
        for i_rct in range(0, n_rct):
            # read reaction
            title.append(f.readline().strip())
            
            t_block = ''
            for i_line in range(0, n_alt_lines):
                t_block = t_block + f.readline()
            t_block = np.array(t_block.split())
            rct[:, i_rct] = t_block.astype(float)
    
    sr1=[]
    sr2=[]
    sp1=[]
    sp2=[]
    sp3=[]
    sp4=[]
    for i_r in range(0, n_rct):
        stmp = title[i_r].split() 
        sr1.append(stmp[0])
        sr2.append(stmp[2])
        sp1.append(stmp[4])
        if (stmp[6] == '+') or (stmp[6] == '-'): # 1 product
            sp2.append('')
            sp3.append('')
            sp4.append('')
        elif (stmp[8] == '+') or (stmp[8] == '-'): # 2 products
            sp2.append(stmp[6])
            sp3.append('')
            sp4.append('')
        elif (len(stmp) <= 10) or (stmp[-1] =='-'): # 3 products
            sp2.append(stmp[6])
            sp3.append(stmp[8])
            sp4.append('')
        else: # 4 products
            sp2.append(stmp[6])
            sp3.append(stmp[8])
            sp4.append(stmp[10])

    reactions = pd.DataFrame({'Reaction':title, 
                              'Reactant1':sr1, 'Reactant2':sr2,
                              'Product1':sp1, 'Product2':sp2, 'Product3':sp3,
                              'Product4':sp4})
    
    rates = {'Altitude':alt, 'Reaction':reactions, 'Rates':rct}
    
    return rates

# Read eflux.out
def eflux(filename):
    alt = []
    with open(filename, 'r') as f:
        # read in number of energy and altitude bins
        n_en, n_alt = map(int, f.readline().split())
        
        # initialize flux array
        flux = np.zeros((n_alt, n_en))
        
        # number of lines in each energy block
        n_lines = int(np.ceil(n_en/10.))
        
        f.readline() # "ENERGY GRID (eV)"
        
        # read in energy grid block
        t_block = ''
        for i_line in range(0, n_lines):
            t_block = t_block + f.readline() 
        t_block = np.array(t_block.split())
        energy = t_block.astype(float)
        
        # read in altitude blocks
        for i_alt in range(0, n_alt):
            # read altitude
            alt.append(f.readline().strip())
            
            t_block = ''
            for i_line in range(0, n_lines):
                t_block = t_block + f.readline()
            # flip negative numbers
            t_block = t_block.replace(' -','  ')
            # there are some missing 'E' if smaller than E-99
            t_block = t_block.replace('-','E-')
            t_block = t_block.replace('EE','E')
            t_block = np.array(t_block.split())
            flux[i_alt, :] = t_block.astype(float)
    
    eflux = {'Altitude':alt, 'Energy':energy, 'Flux':flux}
    
    return eflux