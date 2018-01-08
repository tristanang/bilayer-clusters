import numpy as np
cimport numpy as np

from math import sqrt

def no_z_dist(p1,p2,t):
    p1[2][t] = p2[2][t] = 0

    return np.sum((p1[:,t] - p2[:,t])**2)

def edm_two_loop(matrix1, matrix2=None): #Euclidean distance matrix
	Nconf = matrix1.shape[2]
    Nlipids = matrix1.shape[0]

    if matrix2:
        return
    else:
        dists = np.zeroes([Nlipids,Nlipids,Nconf])

        for t in range(Nconf):
            for i in range(Nlipids):
                for j in range(Nlipids):
                    squared_dist = np.sum
                    dists[i,j,t] =

	
