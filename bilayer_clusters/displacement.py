from bilayer_clusters import constants as c
from bilayer_clusters import boundary as bound

import numpy as np

def linear_gen(start,Nconf,nlog = c.nlog):
    times = []
    t = start
    current = nlog

    t = start + nlog

    while t<Nconf:
        times.append(current)
        current += nlog
        t += nlog

    return np.asarray(times)


def block_displacement(L,lipids):
    return block_displacement_one(L,lipids)

def block_displacement_no(L,lipids): #get displacement of every timestep with reference to the start block
    #L.shape=[Nconf,3]
    v_periodic = np.vectorize(bound.periodic)

    Nconf = lipids.shape[0]
    Nlipids = lipids.shape[1]
    NDIM = lipids.shape[2]

    output_lipids = np.zeros([Nconf,Nlipids,NDIM+1])
    output_lipids[:,:,:3] = lipids

    lipids = lipids[:,:,:2] #lipids.shape is 
    displacement = np.zeros([Nconf,Nlipids,NDIM-1])

    for t in range(Nconf):
        if t%c.nlog == 0:
            frontblock = t
            displacement[frontblock:frontblock+c.nlog,:,:] = lipids[frontblock:frontblock+c.nlog,:,:] - lipids[frontblock,:,:]
        
        else:
            for k in range(NDIM-1):
                displacement[t,:,k] = v_periodic(displacement[t,:,k],L[t,k])

    displacement = displacement**2
    output_lipids[:,:,3] = np.sqrt(displacement[:,:,0]+displacement[:,:,1])

    return output_lipids

def block_displacement_one(L,lipids): #get displacement of every timestep with reference to the start block
    #L.shape=[Nconf,3]
    v_periodic = np.vectorize(bound.periodic)

    Nconf = lipids.shape[0]
    Nlipids = lipids.shape[1]
    NDIM = lipids.shape[2]

    output_lipids = np.zeros([Nconf,Nlipids,NDIM+1])
    output_lipids[:,:,:3] = lipids

    lipids = lipids[:,:,:2]
    
    for t in range(Nconf):
        if t%c.nlog == 0:
            frontblock = t
        else:
            displacement = lipids[t] - lipids[frontblock]
            for k in range(NDIM-1):
                displacement[:,k] = v_periodic(displacement[:,k],L[t,k])

            displacement = displacement**2
            output_lipids[t,:,3] = np.sqrt(displacement[:,0] + displacement[:,1])

    return output_lipids 

def block_displacement_loop(L,lipids):
    Nconf = lipids.shape[0]
    Nlipids = lipids.shape[1]
    NDIM = lipids.shape[2]

    output_lipids = np.zeros([Nconf,Nlipids,NDIM+1])
    output_lipids[:,:,:3] = lipids

    dr = [0,0]

    for t in range(Nconf):
        if t%c.nlog == 0:
            frontblock = t
        else:
            for i in range(Nlipids):
                for k in range(NDIM-1):
                    dr[k] = lipids[t,i,k] - lipids[frontblock,i,k]
                    dr[k] = bound.periodic(dr[k],L[t,k])
                output_lipids[t,i,3] = np.sqrt(dr[0]*dr[0]+dr[1]*dr[1])

    return output_lipids

def s(L,lipids,start,times):
    Nconf = lipids.shape[0]
    Nlipids = lipids.shape[1]
    NDIM = c.NDIM

    v_periodic = np.vectorize(bound.periodic)

    output_lipids = np.zeros([Nconf,Nlipids,NDIM+1])
    output_lipids[:,:,:3] = lipids[:,:,:3]

    if lipids.shape[2] > 3: output_lipids[:,:,3] = lipids[:,:,3]

    for t in times:
        displacement = lipids[t] - lipids[start]
        for k in range(NDIM-1):
            displacement[:,k] = v_periodic(displacement[:,k],L[t,k])

        displacement = displacement**2
        output_lipids[t,:,3] = np.sqrt(displacement[:,0] + displacement[:,1])

    return output_lipids

def linear_displacement(L,lipids,start,Nconf,nlog=c.nlog):
    times = linear_gen(start,Nconf,nlog)
    times += start

    return s(L,lipids,start,times)






