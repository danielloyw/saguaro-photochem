# -*- coding: utf-8 -*-
import numpy as np

#  photo = array of mphotos
#  mphoto = 
#  molname, wav, crs_tot(), crs_ion(), reaction(), brnch(,) 
#  molname  = string
#  wav = 1D array, float
#  crs_tot = 1D array, float
#  crs_ion = 1D array, float
#  brnch = list of (string, 1D array)
#

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
        nbrnch=np.int(snam[1])
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
        
        for z in xrange(0,nbrnch):
            header3 = f.readline()
            sline=''
            for x in xrange(0,8):
                line=f.readline()
                sline=sline+line
    
            sline=np.array(sline.split())
            name.append(snam[0])
            brnch[:,nb]=sline.astype(np.float)
            reaction.append(header3)
            print nb, header3
            nb = nb + 1
        
    f.close()
    
    photo = {"wav":wav,"name":name,"crs_tot":crs_tot,"crs_ion":crs_ion,  \
        "reaction":reaction, "brnch":brnch}

    photo["wav"]=wav
    photo["name"]=name
    photo["crs_tot"]=crs_tot
    photo["crs_ion"]=crs_ion
    photo["reaction"]=reaction
    photo["brnch"]=brnch
    
    return photo   
