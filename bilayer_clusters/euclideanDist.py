import numpy as np
from bilayer_clusters import boundary as bound

from scipy.spatial.distance import cdist as euclidean_distances

"""
def no_z_dist(p1,p2):
    p1[2] = p2[2] = 0

    return np.sum((p1 - p2)**2)

def no_z_squared(p1):
    p1[:,2] = 0
    return np.sum(p1**2,axis=1)

def no_z_xy(p1,p2):
    p1[:,2] = 0 
    p2[:,2] = 0
    return np.dot(p1, p2.T)
"""
def edm_two_loop(L,matrix1, matrix2=None): #Euclidean distance matrix #L array of 3 dimensions
    Nlipids = matrix1.shape[0]

    if not np.any(matrix2):
        matrix2 = matrix1

    Nchol = matrix2.shape[0]

    dists = np.zeros([Nlipids,Nchol])
    r = [0,0]

    for i in range(Nlipids):
        for j in range(Nchol):
            for k in range(2):
                r[k] = matrix1[i][k] - matrix2[j][k]
                r[k] = bound.periodic(r[k],L[k])

            r2 = np.sqrt(r[0]*r[0] + r[1]*r[1])
            dists[i,j] = r2
    
    return dists

def edm(L,matrix1, matrix2=None): #Euclidean distance matrix
    Nlipids = matrix1.shape[0]

    if not np.any(matrix2):
        matrix2 = matrix1

    Nchol = matrix2.shape[0]

    NDIM = 3
    lst = [0,0]
    dists = np.zeros([Nlipids,Nchol])

    #v_periodic = bound.periodic
    v_periodic = np.vectorize(bound.periodic)

    for k in range(NDIM-1):
        matrix_1 = matrix1[:,k].reshape(-1,1)
        matrix_2 = matrix2[:,k].reshape(-1,1)
        dist1d = euclidean_distances(matrix_1,matrix_2,'euclidean')
        dist1d = v_periodic(dist1d,L[k])
        lst[k] = dist1d

    dists = np.sqrt(lst[0]**2 + lst[1]**2)

    #dists = euclidean_distances(matrix1,matrix2,'euclidean')

    """
    x_squared = no_z_squared(matrix1[t])
    y_squared = no_z_squared(matrix2[t])
    xy = no_z_xy(matrix1[t],matrix2[t])
    dists[t] = np.sqrt(-2*xy + y_squared + x_squared)
    """
    
    return dists
    
