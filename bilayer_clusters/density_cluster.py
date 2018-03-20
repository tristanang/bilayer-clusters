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

def randomCluster(siz,arr):
    Nparticles = len(arr)
    idx = np.random.randint(Nparticles, size=siz)

    cluster = arr[idx,:]

    return cluster