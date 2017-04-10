# -*- coding: utf-8 -*-
import numpy as np

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
            reaction.append(header3)
            print nb, header3
            nb = nb + 1
        
    f.close()
    
    photo = {"name":name,"wav":wav,"crs_tot":crs_tot,"crs_ion":crs_ion,  \
        "reaction":reaction, "brnch":brnch}

    photo["name"]=name
    photo["wav"]=wav
    photo["crs_tot"]=crs_tot
    photo["crs_ion"]=crs_ion
    photo["reaction"]=reaction
    photo["brnch"]=brnch
    
    return photo   