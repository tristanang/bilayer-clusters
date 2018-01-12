import numpy as np 
from bilayer_clusters.gvf import *

def test_gvf1():
    cluster = [np.array([1,1,1,1]),np.array([5,5,5,5,5]),np.array([9,9,9,9])]
    assert gvf_helper(cluster) == 1

    cluster = [np.arange(1,4),np.arange(100,104),np.arange(10000,10004)]

    print(gvf_helper(cluster))

    return True


if __name__ == '__main__':
    test_gvf1()
    print("passed")