# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 11:56:25 2014

@author: rogeryelle
"""

# -*- coding: utf-8 -*-
import numpy as np
from csv_to_list import csv_to_list

# script to read atmosphere files

def read_atm(filename):
    
    f=open(filename,'r')
        
    nalt, nmol = map(int,f.readline().split())
    cos_sza = np.float(f.readline().strip())

    line=f.readline()    
    sline='Ntot'
    for x in xrange(0,32):
        sline=sline+f.readline()       
    name=sline.split()
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        line=f.readline()
        sline=sline+line 
    sline=np.array(sline.split())
    alt=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline()
    sline=np.array(sline.split())
    rad=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    grv=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    Tn=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    Te=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    prs=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    rho=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    mmw=sline.astype(np.float)
    
    line=f.readline()
    sline=''
    for x in xrange(0,31):
        sline=sline+f.readline() 
    sline=np.array(sline.split())
    Kzz=sline.astype(np.float)
      
    den=np.zeros((nalt,nmol+1), dtype=np.float)
    for n in xrange(0,nmol+1):
        line=f.readline()
        sline=''
        for x in xrange(0,31):
            sline=sline+f.readline() 
        sline=np.array(sline.split())
        den[:,n]=sline.astype(np.float)

    atm = {"nalt":nalt, "nmol":nmol, "cos_sza":cos_sza,"name":name,   \
        "alt":alt, "rad":rad, "grv":grv, "Tn":Tn, "Te":Te, "prs":prs, \
        "rho":rho, "mmw":mmw, "Kzz":Kzz, "den":den}
    atm["nalt"] = nalt
    atm["nmol"] = nmol
    atm["cos_sza"] = cos_sza
    atm["name"] = name
    atm["alt"] = alt
    atm["rad"] = rad
    atm["Tn"] = Tn
    atm["Te"] = Te
    atm["prs"] = prs
    atm["rho"] = rho
    atm["mmw"] = mmw
    atm["Kzz"] = Kzz
    atm["den"] = den
   
    f.close()
    return atm

# script to read photoA.dat

def read_photoA(filename):
    
    f=open(filename,'r')
    
    header1 = f.readline()
    header1 = header1.strip()
    [nwav,nmol,nchn]=header1.split()
    nwav=np.int(nwav)
    nmol=np.int(nmol)
    nchn=np.int(nchn)
    print nwav,nmol,nchn
    crs_tot = np.zeros([nwav,nmol])
    crs_ion = np.zeros([nwav,nmol])
    name=[]
    reaction=[]
    nbrnch=np.zeros([nmol], dtype='int8')
    brnch = np.zeros([nwav,nchn])
    header2 = f.readline()
    header2 = header2.strip()
    
#    pmol = np.zeros((2,),dtype=('a10,f4,a10'))    
    
    sline=''
    for x in xrange(0,8):
        line=f.readline()
        sline=sline+line
    
    sline=np.array(sline.split())
    wav=sline.astype(np.float)
    
    nb = 0
    for y in xrange(0,nmol):
        header1 = f.readline()
        header1 = header1.strip()
        snam = header1.split()
        name.append(snam[0])        
        nbrnch[y]=np.int(snam[1])
        header2 =f.readline()
        
        sline=''
        for x in xrange(0,8):
            line=f.readline()
            sline=sline+line
    
        sline=np.array(sline.split())
        crs_tot[:,y]=sline.astype(np.float)
        
        header2 = f.readline()
        print header2
    
        sline=''
        for x in xrange(0,8):
            line=f.readline()
            sline=sline+line
    
        sline=np.array(sline.split())
        crs_ion[:,y]=sline.astype(np.float)
        
        for z in xrange(0,nbrnch[y]):
            header3 = f.readline()
            sline=''
            for x in xrange(0,8):
                line=f.readline()
                sline=sline+line
    
            sline=np.array(sline.split())
            brnch[:,nb]=sline.astype(np.float)
            reaction.append(header3)
            print nb, header3
            nb = nb + 1
        
    f.close()
    
#    photo = {"wav":wav,"name":name,"crs_tot":crs_tot,"crs_ion":crs_ion,  \
#        "reaction":reaction, "brnch":brnch}

#    photo["wav"]=wav
#    photo["name"]=name
#    photo["crs_tot"]=crs_tot
#    photo["crs_ion"]=crs_ion
#    photo["reaction"]=reaction
#    photo["brnch"]=brnch
    
#    return photo   
    return wav, name, crs_tot, crs_ion, nbrnch, reaction, brnch

# script to read photoC.dat 
   
def read_photoC(filename):
    
    f=open(filename,'r')
    
    header1 = f.readline()
    header1 = header1.strip()
    [nwav,nmol,nchn]=header1.split()
    nwav=np.int(nwav)
    nmol=np.int(nmol)
    nchn=np.int(nchn)
    print nwav,nmol,nchn
    crs_tot = np.zeros([nwav,nmol])
    crs_ion = np.zeros([nwav,nmol])
    name=[]
    reaction=[]
    brnch = np.zeros([nwav,nchn])
    header2 = f.readline()
    header2 = header2.strip()
    
    sline=''
    for x in xrange(0,20):
        line=f.readline()
        sline=sline+line
    
    sline=np.array(sline.split())
    wav=sline.astype(np.float)
    
    nb = 0
    for y in xrange(0,nmol):
        header1 = f.readline()
        header1 = header1.strip()
        snam = header1.split()
        nbrnch = np.int(snam[1])
        print name, nbrnch
        header2 =f.readline()
        print header2
        
        sline=''
        for x in xrange(0,20):
            line=f.readline()
            sline=sline+line
    
        sline=np.array(sline.split())
        crs_tot[:,y]=sline.astype(np.float)
        
        header2 = f.readline()
        print header2
    
        sline=''
        for x in xrange(0,20):
            line=f.readline()
            sline=sline+line
    
        sline=np.array(sline.split())
        crs_ion[:,y]=sline.astype(np.float)
        
        for z in xrange(0,nbrnch):
            header3 = f.readline()
            sline=''
            for x in xrange(0,20):
                line=f.readline()
                sline=sline+line
    
            sline=np.array(sline.split())
            brnch[:,nb]=sline.astype(np.float)
            name.append(snam[0])
            reaction.append(header3.strip())
            print nb, header3
            nb = nb + 1
        
    f.close()
    
    photo = {"wav":wav,"name":name,"crs_tot":crs_tot,"crs_ion":crs_ion,  \
        "nbrnch":nbrnch,"reaction":reaction, "brnch":brnch}

    photo["wav"]=wav
    photo["name"]=name
    photo["crs_tot"]=crs_tot
    photo["crs_ion"]=crs_ion
    photo["nbrnch"]=nbrnch
    photo["reaction"]=reaction
    photo["brnch"]=brnch
    
    return photo       
    
# script to read nmolecules.dat
    
def read_nmolecules(filename):
    
    f=open(filename,'r')
    
    header1 = f.readline()
    nmol=np.int(header1.strip())
    slin = f.readline()
    name=[]
    stat=np.zeros(nmol, dtype=np.int)
    chrg=np.zeros(nmol, dtype=np.int)
    wght=np.zeros(nmol)
    hyd=np.zeros(nmol, dtype=np.int)
    car=np.zeros(nmol, dtype=np.int)
    n14=np.zeros(nmol, dtype=np.int)
    n15=np.zeros(nmol, dtype=np.int)
    oxy=np.zeros(nmol, dtype=np.int)
    for x in xrange(0,136):
        slin=f.readline()
        slin=slin.split()
        name.append(slin[1])
        stat[x]=np.int(slin[2])
        chrg[x]=np.int(slin[3])
        wght[x]=np.float(slin[4])
        hyd[x]=np.int(slin[5])
        car[x]=np.int(slin[6])
        n14[x]=np.int(slin[7])
        n15[x]=np.int(slin[8])
        oxy[x]=np.int(slin[9])
        
    f.close()
    
    molec={"name":name,"stat":stat,"chrg":chrg,"wght":wght,"hyd":hyd, \
        "car":car,"n14":n14,"n15":n15,"oxy":oxy}
    molec["name"]=name
    molec["stat"]=stat
    molec["chrg"]=chrg
    molec["wght"]=wght
    molec["hyd"]=hyd
    molec["car"]=car
    molec["n14"]=n14
    molec["n15"]=n15
    molec["oxy"]=oxy
    
    return molec          

# scipt to read nreactions.csv 
    
def read_reactions(filename):

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

    return react   