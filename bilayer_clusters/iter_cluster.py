import numpy as np
import random
from functools import reduce

from bilayer_clusters import euclideanDist
from bilayer_clusters import constants as c
from bilayer_clusters import density_cluster as dc

def combine(lipids,chol):
    Nlipids = len(lipids)
    Nchol = len(chol)
    Ntotal = Nlipids+Nchol

    lst1 = np.zeros([Nlipids,1])
    lipids = np.concatenate((lipids,lst1),axis=1)

    lst1 = np.zeros([Nchol,1])
    lst1 += 1
    chol = np.concatenate((chol,lst1),axis=1)

    return np.concatenate((lipids,chol),axis=0)

def cluster_labels(flag,L,total):
    Ntotal = len(total)

    identity_arr = total[:,-1]
    identity_arr = identity_arr.astype(int)

    labels = np.zeros(Ntotal,dtype=int)
    labels -= 1

    neighbors = [[] for i in range(Ntotal)]


    cutoff = c.cutoff[flag]

    edm = euclideanDist.edm(L,total)

    for i in range(Ntotal):
        for j in range(i+1,Ntotal):
            if edm[i,j] <= cutoff[identity_arr[i]+identity_arr[j]]:
                neighbors[i].append(j)


    clusters = []

    for i in range(Ntotal):
        group = set_neighbors(i,neighbors)
        group.add(i)

        if clusters == []: 
            clusters.append(group)
        else:
            i = 0
            while i < len(clusters):
                cluster = clusters[i]
                if group.intersection(cluster):
                    group |= clusters.pop(i)
                else:
                    i += 1
            clusters.append(group)

    nclust = 0

    for group in clusters:
        if len(group) != 1:

            for elem in group:
                labels[elem] = nclust

            nclust += 1
    """
    for i in range(Ntotal):
        if labels[i] == -1:
            labels[i] = nclust
            cluster = set_neighbors(i,neighbors)
            for elem in cluster:
                labels[elem] = nclust

            #set neighbors

            nclust += 1
    """
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

cluster_sizes = dc.cluster_sizes
means = dc.means

def randomCluster(Nlipids,Nchol,lipids,chol):
    Nparticles = len(arr)
    idx = np.asarray(random.sample(list(range(Nparticles)),siz))

    cluster = arr[idx,:]

    return cluster

def mean_cluster_size(total,L,flag='upper'):
    labels = cluster_labels(flag,L,total)
    hsh = cluster_sizes(labels)
    return means(hsh)

def meanRandom(originalArr,L,size,flag='upper',iter=10):
    lst = []

    for i in range(iter):
        cluster = randomCluster(size,originalArr)
        lst.append(mean_cluster_size(cluster,L,flag))

    lst=list(zip(*lst))

    mean=sum(lst[0])/iter
    mS=sum(lst[1])/iter
    
    return mean,mS

def normSize(total,L,originalArr,flag='upper',iter=10,):
    nmean,nmS = mean_cluster_size(total,L,flag)
    #print(nmean,nmS)
    #denominator
    size = len(arr)
    dmean,dmS = meanRandom(originalArr,L,size,flag,iter)

    return (nmean/dmean,nmS/dmS)

#def counter(cluster):
    


