
def read_nmolecules(filename):
    import numpy as np    
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