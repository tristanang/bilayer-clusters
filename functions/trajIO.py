import numpy as np
cimport numpy as np

import boundary as bound

import re

def cholConc(filename):
    file = open(filename,"r")
    for line in file:
        m = re.match("CHOL", line)
        if m:
            line = line.split()
            file.close()
            return int(line[1]) * 2

cdef processTrajCOM(str trajFileName,int Nchol,int NDIM,int Nconf):
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
    
    L = np.zeros((NDIM, Nconf), dtype=np.double)
    com_lipids = np.zeros((Nlipids,NDIM,Nconf), dtype=np.double)
    com_chol = np.zeros((Nchol,NDIM,Nconf), dtype=np.double)
    
    x = np.zeros((Nlipidbeads,NDIM), dtype=np.double)
    y = np.zeros((Ncholbeads,NDIM), dtype=np.double)
    
    cdef int t,i,j,k,index
    
    for t in range(Nconf):
        trajFile.readline()
        
        #Box Sizes
        for k in range(NDIM):
            L[k,t] = float(trajFile.readline().strip())
        
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
                    x[ii+j,k] = bound.periodic_p(x[ii+j,k],x[ii,k],L[k,t])
                    com_lipids[i,k,t] += x[ii+j,k]/Nperlipid
            
            for i in range(Nchol):
                ii = i * Nperchol
                for j in range(Nperchol):
                    y[ii+j,k] = bound.periodic_p(y[ii+j,k],y[ii,k],L[k,t])
                    com_chol[i,k,t] += y[ii+j,k]/Nperchol
    
    trajFile.close()

    return N,L,com_lipids,com_chol