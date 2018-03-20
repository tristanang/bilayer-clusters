import numpy as np 
from bilayer_clusters import euclideanDist
from functools import reduce
import matplotlib.pyplot as plt

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
        self.L = L
        self.cutoff = cutoff
        #self.arr = arr
        self.Nparticles = len(arr)
        edm = euclideanDist.edm(L,arr)
        #edm = euclidean_distances(arr,arr,'euclidean')
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

        self.sizes = self.getSizes()

    """   
    def singletonMerge(self):
        singletons = list(filter(lambda x : len(x) == 1,self.clusters))
        self.clusters = list(filter(lambda x : len(x) != 1,self.clusters))
        
        if singletons:
            merge = reduce(lambda x,y : x | y,singletons)
            self.clusters.append(merge)
        return
    """

    def getCenter(self):
        return

    def normalize(self,original,trials=5):
        def hidden_r(x):
            random_clust = randomCluster(self.Nparticles,original)
            c = Cluster(random_clust,self.L,self.cutoff)
            return c.sizes 

        temp = [0 for i in range(trials)]
        temp = list(map(hidden_r, temp))

        m = list(map(np.mean,temp))

        temp = list(map(np.square,temp))
        mS = list(map(np.mean,temp))

        """
        #parallelize
        for i in range(trials):
            random_clust = randomCluster(self.Nparticles,original)
            c = Cluster(random_clust,self.L,self.cutoff)
            sizes = c.sizes
            
            m[i] = np.mean(sizes)

            sizes = sizes ** 2
            mS[i] = np.mean(sizes)
        """

        return (self.mean()/np.mean(m)),(self.weightedMean()/(np.mean(mS)/np.mean(m)))



    def getSizes(self):
        #output = list(map(lambda x: float(len(x))/self.Nparticles, self.clusters))
        output = pool.map(len, self.clusters)
        output = np.asarray(output)

        return output

    """
    def hidden_f(self,sett):
        output = []
        for i in sett:
            output.append((self.arr[i,0],self.arr[i,1]))

        return output

    def spitXY(self): ##
        return list(map(self.hidden_f,self.clusters))
    """

    def test(self):
        assert (sum(self.getSizes())-self.Nparticles) <= 0.0000001

    def weightedMean(self):
        arr = self.sizes ** 2
        return np.mean(arr)/(self.mean())

    def mean(self):
        return np.mean(self.sizes)

def randomCluster(siz,arr):
    Nparticles = len(arr)
    idx = np.random.randint(Nparticles, size=siz)

    cluster = arr[idx,:]

    return cluster