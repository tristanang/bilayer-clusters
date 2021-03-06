import numpy as np
cimport numpy as np
cimport cython

from bilayer_clusters import boundary as bound
import constants as c

import re

def cholConc(filename):
    file = open(filename,"r")
    for line in file:
        m = re.match("CHOL", line)
        if m:
            line = line.split()
            file.close()
            return int(line[1]) * 2

def processTrajCOM(str trajFileName,int Nchol,int NDIM,int Nconf):
    trajFile = open(trajFileName,'r')
    cdef int N = int(trajFile.readline().split()[0])
    cdef int Nperlipid = c.Nperlipid
    cdef int Nperchol = c.Nperchol
    cdef int Nlipids = (N - Nperchol * Nchol) // Nperlipid
    cdef int Ncholbeads = Nchol * Nperchol
    cdef int Nlipidbeads = Nlipids * Nperlipid
    
    #reading traj file
    trajFile.close()
    trajFile = open(trajFileName,'r')

    cdef np.ndarray[np.double_t, ndim=2] L,x,y
    cdef np.ndarray[np.double_t, ndim=3] com_lipids, com_chol
    
    L = np.zeros((Nconf,NDIM), dtype=np.double)
    com_lipids = np.zeros((Nconf,Nlipids,NDIM), dtype=np.double)
    com_chol = np.zeros((Nconf,Nchol,NDIM), dtype=np.double)
    
    x = np.zeros((Nlipidbeads,NDIM), dtype=np.double)
    y = np.zeros((Ncholbeads,NDIM), dtype=np.double)
    
    cdef int t,i,j,k,index
    
    for t in range(Nconf):
        trajFile.readline()
        
        #Box Sizes
        for k in range(NDIM):
            L[t,k] = float(trajFile.readline().strip())
        
        #Bead Coordinates
        for i in range(Nlipidbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
                
            for k in range(NDIM):
                x[i,k] = float(config[k])
            
        for i in range(Ncholbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
            
            for k in range(NDIM):
                y[i,k] = float(config[k])
            
        for i in range(Nlipidbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
             
            for k in range(NDIM):
                index = i+Nlipidbeads//2
                x[index,k] = float(config[k])
        
        for i in range(Ncholbeads//2):
            line = trajFile.readline().split()
            config = [line[j] for j in range(1,NDIM+1)]
            
            for k in range(NDIM):
                index = i+Ncholbeads//2
                y[index,k] = float(config[k])
    
        for k in range(NDIM):
            for i in range(Nlipids):
                ii = i * Nperlipid
                for j in range(Nperlipid):
                    x[ii+j,k] = bound.periodic_p(x[ii+j,k],x[ii,k],L[t,k])
                    com_lipids[t,i,k] += x[ii+j,k]/Nperlipid
            
            for i in range(Nchol):
                ii = i * Nperchol
                for j in range(Nperchol):
                    y[ii+j,k] = bound.periodic_p(y[ii+j,k],y[ii,k],L[t,k])
                    com_chol[t,i,k] += y[ii+j,k]/Nperchol
    
    trajFile.close()

    return N,L,com_lipids,com_chol

def translateZ(np.ndarray[np.double_t, ndim=3] x,np.ndarray[np.double_t, ndim=3] y):
    cdef int Nconf,Nlipidbeads,Ncholbeads,N
    cdef double z_avg
    cdef int t,i
    
    Nconf = x.shape[0]
    Nlipidbeads = x.shape[1]
    Ncholbeads = y.shape[1]
    N = Nlipidbeads + Ncholbeads

    for t in range(Nconf):
        z_avg = 0.
            
        for i in range(Nlipidbeads):
            z_avg += x[t,i,2]

        for i in range(Ncholbeads):
            z_avg += y[t,i,2]

        z_avg = z_avg / N

        for i in range(Nlipidbeads):
            x[t,i,2] -= z_avg

        for i in range(Ncholbeads):
            y[t,i,2] -= z_avg

    return x,y

def layering(com_lipids):
    upper_lipids = com_lipids[com_lipids[:,2]>=0] 
    lower_lipids = com_lipids[com_lipids[:,2]<0]

    return upper_lipids,lower_lipids

def decompress(file="comTraj.npz"):
    npzfile = np.load(file)

    L = npzfile['L']
    com_lipids = npzfile['com_lipids']
    com_chol = npzfile['com_chol']

    return L,com_lipids,com_chol

def rawOrCOM(filename):
    return False


