import numpy as np
import random
from functools import reduce

from bilayer_clusters import euclideanDist
from bilayer_clusters import constants as c

def cluster_labels(flag,L,lipids,chol=[]):
    Nlipids = len(lipids)
    Nchol = len(chol)
    Ntotal = Nlipids+Nchol

    identity_arr = np.zeros(Nlipids+Nchol,dtype=int)
    identity_arr[Nlipids:] = 1

    labels = np.zeros(Nlipids+Nchol,dtype=int)
    labels -= 1

    neighbors = [[] for i in range(Ntotal)]

    assert np.count_nonzero(identity_arr) == Nchol

    cutoff = c.cutoff[flag]

    if Nchol:
        total = np.concatenate((lipids, chol), axis=0)
    else:
        total = chol

    edm = euclideanDist.edm(L,total)

    for i in range(Ntotal):
        for j in range(i+1,Ntotal):
            if edm[i,j] <= cutoff[identity_arr[i]+identity_arr[j]]:
                neighbors[i].append(j)

    nclust = 0

    for i in range(Ntotal):
        if labels[i] == -1:
            labels[i] = nclust
            cluster = set_neighbors(i,neighbors)
            for elem in cluster:
                labels[elem] = nclust

            #set neighbors

            nclust += 1

    return labels

def set_neighbors(i,neighbors):
    
    assert type(i) == type(1)

    cluster = set(neighbors[i])
    #print(cluster)
    if neighbors[i] == []:
        return cluster
    else:
        for elem in neighbors[i]:
            cluster |= set_neighbors(elem,neighbors)

    return cluster
        #cluster | reduce((lambda x,y: set_neighbors(x,neighbors) | set_neighbors(y,neighbors)),cluster)


