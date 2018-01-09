import sys
sys.path.append("../")

from bilayer_clusters import trajIO
from bilayer_clusters.euclideanDist import *

import numpy as np
import random
"""
def test_no_z_dist():
    p1 = np.random.rand(3)
    p2 = np.random.rand(3)

    dist1 = 0

    for k in range(2):
        dist1 += (p1[k] - p2[k])**2

    dist2 = no_z_dist(p1,p2)

    if dist1 == dist2:
        print("no_z_dist function passed")
        return True
    else:
        print("no_z_dist function failed")
        return False
"""
def test_2vs0_loop():
    N = 23
    N2 = 43
    NDIM = 3
    L = [1,1,1]

    p1 = np.random.rand(N,NDIM)
    p2 = np.random.rand(N2,NDIM)

    dist1 = edm_two_loop(L,p1,p2)
    dist2 = edm(L,p1,p2)

    sum = 0.0

    for t in range(1):
        difference = np.linalg.norm(dist1 - dist2, ord='fro')
        sum += difference

    if sum < 0.00000001:
        print(str(sum)+" 2loop vs 0 loop passed")
        return True
    else:
        print("2loop vs 0 loop failed")
        return False

if __name__ == "__main__":
    test_2vs0_loop()
