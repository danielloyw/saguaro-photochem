
def read_molecules(filename):
    import numpy as np    
    f=open(filename,'r')
    
    f.readline()
    f.readline()
    name=[]
    stat=[]
    chrg=[]
    wght=[]
    hyd=[]
    car=[]
    n14=[]
    n15=[]
    oxy=[]

    for line in f:
        line=line.split()
        name.append(line[1])
        stat.append(np.int(line[2]))
        chrg.append(np.int(line[3]))
        wght.append(np.float(line[4]))
        hyd.append(np.int(line[5]))
        car.append(np.int(line[6]))
        n14.append(np.int(line[7]))
        n15.append(np.int(line[8]))
        oxy.append(np.int(line[9]))
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