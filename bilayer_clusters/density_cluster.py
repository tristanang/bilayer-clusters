from sklearn.cluster import DBSCAN
from scipy.cluster import hierarchy

import numpy as np

from bilayer_clusters import euclideanDist 

def dbscan_wrapper(arr,L,cutoff):
    edm = euclideanDist.edm(L,arr)
    db = DBSCAN(eps=cutoff, min_samples=2, metric='precomputed', n_jobs=-1)
    #db.fit(edm)
    #labels = db.labels_
    labels = db.fit_predict(edm)

    return labels

def hierarchy_wrapper(arr,L,cutoff):
    edm = euclideanDist.edm(L,arr)
    y = edm[np.triu_indices(edm.shape[0],1)]
    link = hierarchy.linkage(y, method='single')
    cut = hierarchy.cut_tree(link,height=cutoff)
    cut = [y for x in cut for y in x]

    return np.asarray(cut)

def cluster_sizes(labels):
    unique, counts = np.unique(labels, return_counts=True)
    return dict(zip(unique, counts))

def means(hsh):
    singletons = hsh.pop(-1,0.0)
    number_clusters = singletons + len(hsh)

    lst = np.asarray(list(hsh.values()))

    sum_squared = singletons + sum(lst**2) 
    sum_cluster = sum(lst) + singletons

    mean = float(sum_cluster)/number_clusters
    mS = float(sum_squared)/sum_cluster

    return mean,mS

def mean_cluster_size(arr,L,cutoff,alg=dbscan_wrapper):
    labels = alg(arr,L,cutoff)
    hsh = cluster_sizes(labels)
    return means(hsh)

def meanRandom(originalArr,L,cutoff,size,iter=10 ,alg=dbscan_wrapper):
    lst = []

    for i in range(iter):
        cluster = randomCluster(size,originalArr)
        lst.append(mean_cluster_size(cluster,L,cutoff))

    lst=list(zip(*lst))

    mean=sum(lst[0])/iter
    mS=sum(lst[1])/iter

    return mean,mS

def normSize(arr,L,cutoff,originalArr,alg=dbscan_wrapper):
    #numerator
    nmean,nmS = mean_cluster_size(arr,L,cutoff)
    #print(nmean,nmS)
    #denominator
    size = len(arr)
    dmean,dmS = meanRandom(originalArr,L,cutoff,size)
    print(dmean,dmS)

    return (nmean/dmean,nmS/dmS)

def randomCluster(siz,arr):
    Nparticles = len(arr)
    idx = np.random.randint(Nparticles, size=siz)

    cluster = arr[idx,:]

    return cluster