# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 14:56:25 2014

@author: jterwin
"""

from pylab import *



def read_cks(fname):
    """reads in the correlated k coefficient in file fname
    returns a dictionary
    
    ex.
    >>> cks = read_cks('co2.ckl')"""
    f = open(fname)
    ng, nwav, ntmp, nprs = map(int,f.readline().split())
    
    w = array(map(float,f.readline().split()))

    tmp = zeros(ntmp)
    prs = zeros(nprs)
    wav = zeros(nwav+1)
    coeffs = zeros((ng,nwav,ntmp,nprs))

    # loop to read in tmp's prs's, wav's, and coeffs
    for nt in range(ntmp):
        for np in range(nprs):
            #print 't,p = ', f.readline().split()
            tmp[nt],prs[np] = map(float,f.readline().split())
            
            for nw in range(nwav):
                #print 'info = ', f.readline().split()
                vals = map(float,f.readline().split())
                wav[nw] = vals[0]
                wav[nw+1] = vals[1]
                coeffs[:,nw,nt,np] = vals[2:2+ng]
    
    # put in dictionary to return
    cks = {"ng":ng, "nwav":nwav, "ntmp":ntmp, "nprs":nprs}
    cks["w"] = w
    cks["tmp"] = tmp
    cks["prs"] = prs
    cks["wav"] = wav
    cks["coeffs"] = coeffs
    cks["xck"] = array([0.032076, 0.157791, 0.341886, 0.525981, 0.651696, 0.693915, 0.733670, 0.791886, 0.850102, 0.889857, 0.903208, 0.915779, 0.934189, 0.952598, 0.965170, 0.969392, 0.973367, 0.979189, 0.985010, 0.988986, 0.990321, 0.991578, 0.993419, 0.995260, 0.996517, 0.996986, 0.997567, 0.998419, 0.999270, 0.999852])
    
    return cks
    

def sigma_vs_nu_plot(cks,ax=None):
    """ plots the expected value of sigma vs nu
    if you don't provide ax, then current axis is used
    
    future: be able to select P,T or np,nt
    """
    from matplotlib import pyplot as plt
    import numpy
    
    nt = 0
    np = 0
    print cks["prs"][np], cks["tmp"][nt]    
    
    wav = (cks["wav"][:-1]+cks["wav"][1:])/2.0
    sig = numpy.zeros(cks["nwav"])
    for nw in range(cks["nwav"]):
        sig[nw] = cks["w"].dot(cks["coeffs"][:,nw,np,nt])
    
    if ax is None:
        ax = plt.gca()
        
    return ax.semilogy(wav,sig)


'''
cks1 = read_cks('co2_temp1.ckl')
cks2 = read_cks('co2_temp2.ckl')

fig, ax = subplots()
sigma_vs_nu_plot(cks1,ax)
sigma_vs_nu_plot(cks2,ax)
'''


cks2 =  read_cks('data/co2_3450_3750_am.ckl')
#cks4 = read_cks('data/co2_4_ry.ckl')
cks4 =  read_cks('data/co2_2200_2400_am.ckl')
cks15 = read_cks('data/co2_500_850_am.ckl')

#for key in cks15.keys():
#    print key, cks15[key]


fig, ax = subplots()
sigma_vs_nu_plot(cks2)
sigma_vs_nu_plot(cks4)
sigma_vs_nu_plot(cks15)

fig, axs = subplots(1,3,sharey=True)
sigma_vs_nu_plot(cks2,axs[2])
sigma_vs_nu_plot(cks4,axs[1])
sigma_vs_nu_plot(cks15,axs[0])




'''

# plot 
cks = cks1
nw = 5
nt = abs(cks["tmp"]-150.0).argmin()
np = abs(cks["prs"]-1.0e1).argmin()

for nw in range(cks["nwav"]):
    figure(100+nw)

    for np in range(cks["nprs"]):
        if any(cks["coeffs"][:,nw,nt,np]>0):
            semilogy(cks["xck"],cks["coeffs"][:,nw,nt,np], label=str(np))

    legend(loc=2)
    title(str(nw))

cks = cks2
nw = 5
nt = abs(cks["tmp"]-150.0).argmin()
np = abs(cks["prs"]-1.0e1).argmin()

for nw in range(cks["nwav"]):
    figure(100+nw)

    for np in range(cks["nprs"]):
        if any(cks["coeffs"][:,nw,nt,np]>0):
            semilogy(cks["xck"],cks["coeffs"][:,nw,nt,np], label=str(np))

    legend(loc=2)
    title(str(nw))

'''

show()
