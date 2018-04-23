from bilayer_clusters import iter_cluster
from bilayer_clusters import trajIO
import numpy as np

def test_findCluster():
    lst = []
    set1 = set([i for i in range(10)])
    set2 = set([i for i in range(14,29)])
    lst = [set1,set2]

    assert iter_cluster.findCluster(lst,5) == 0
    assert iter_cluster.findCluster(lst,20) == 1
    assert iter_cluster.findCluster(lst,1000) == None

"""
def test_buildCluster():
    lst = [(i,i) for i in range(100)]
    lst.append((0.1,0.1))

    c = iter_cluster.Cluster(lst,[],0.5)
    #print(c.clusters)
    c.singletonMerge()
    print(c.clusters)
    print("-----")
    print(c.Nparticles)
    print(c.getSizes())
"""

def test_buildCluster(): #& singletonMerge
    L, com_lipids, com_chol = trajIO.decompress("comTraj.npz")
    com_lipids, com_chol = trajIO.translateZ(com_lipids, com_chol)
    Nlipids = len(com_lipids[0])

    c = iter_cluster.Cluster(com_lipids[0],L[0],0.3)
    c.test()
    """
    manualSingleton = set()

    for cluster in c.clusters:
        if len(cluster) == 1:
            manualSingleton |= cluster

    c.singletonMerge()
    c.test()

    assert c.clusters[-1] == manualSingleton
    """
    return

def test_randomCluster():
    L, com_lipids, com_chol = trajIO.decompress("comTraj.npz")
    com_lipids, com_chol = trajIO.translateZ(com_lipids, com_chol)
    arr = com_lipids[0]
    size = 43

    cluster = iter_cluster.randomCluster(size,arr)
    assert len(cluster) == 43


    return

def test_randomIterCluster():
    L, com_lipids, com_chol = trajIO.decompress("comTraj.npz")
    com_lipids, com_chol = trajIO.translateZ(com_lipids, com_chol)
    size = 45


    return

if __name__ == '__main__':
    test_findCluster()
    test_buildCluster()
    test_randomCluster()
    test_randomIterCluster()
    print("Passed!")
