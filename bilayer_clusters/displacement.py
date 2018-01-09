from bilayer_clusters import constants as c
from bilayer_clusters import boundary as bound

import numpy as np

def block_displacement(L,lipids): #get displacement of every timestep with reference to the start block
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
            displacement = lipids[frontblock:frontblock+c.nlog,:,:] - lipids[frontblock,:,:]
            
            for dt in range(c.nlog):
                for k in range(NDIM-1):
                    displacement[dt,:,k] = v_periodic(displacement[dt,:,k],L[t,k])

            displacement = displacement**2
            output_lipids[frontblock:frontblock+c.nlog,:,3] = np.sqrt(displacement[:,:,0] + displacement[:,:,1])

        else:
            continue    

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




