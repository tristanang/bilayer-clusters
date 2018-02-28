import numpy as np 
from bilayer_clusters import euclideanDist
from functools import reduce
import matplotlib.pyplot as plt

from scipy.spatial.distance import cdist as euclidean_distances

def findCluster(clusters,item):
    Nclusters = len(clusters)

    for i in range(Nclusters):
        if item in clusters[i]:
            return i

    return None

def display(cluster):
    colors = []

    for group in cluster:
        x, y = zip(*group)
        plt.plot(x, y)
    return

class Cluster:
    def __init__(self,arr,L,cutoff):
        self.arr = arr
        self.Nparticles = len(arr)
        #edm = euclideanDist.edm(L,arr)
        edm = euclidean_distances(arr,arr,'euclidean')
        self.clusters = [set([i]) for i in range(self.Nparticles)]

        for i in range(self.Nparticles):
            cut = edm[i,:]

            indices = np.ndarray.flatten(np.argwhere(cut <= cutoff))
            indices = set(map(lambda x : findCluster(self.clusters, x),indices))
            indices = sorted(list(indices), reverse = True)
            small = indices.pop()

            for index in indices:
                self.clusters[small] |= self.clusters[index]
                del self.clusters[index]
          
    def singletonMerge(self):
        singletons = list(filter(lambda x : len(x) == 1,self.clusters))
        self.clusters = list(filter(lambda x : len(x) != 1,self.clusters))
        
        if singletons:
            merge = reduce(lambda x,y : x | y,singletons)
            self.clusters.append(merge)
        return

    def getCenter(self):
        return

    def getSizes(self):
        output = list(map(lambda x: float(len(x))/self.Nparticles, self.clusters))

        return output

    def hidden_f(self,sett):
        output = []
        for i in sett:
            output.append((self.arr[i,0],self.arr[i,1]))

        return output

    def spitXY(self): ##
        return list(map(self.hidden_f,self.clusters))

    def test(self):
        assert (sum(self.getSizes())-1) <= 0.0000001

def randomCluster(siz,arr):
    Nparticles = len(arr)
    idx = np.random.randint(Nparticles, size=siz)

    cluster = arr[idx,:]

    return cluster