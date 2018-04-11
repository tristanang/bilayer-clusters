import sys
sys.path.append("../")

from bilayer_clusters import trajIO
from bilayer_clusters.euclideanDist import *

import numpy as np
import random

def time_function(f, *args):
  """
  Call a function f with args and return the time (in seconds) that it took to execute.
  """
  import time
  tic = time.time()
  f(*args)
  toc = time.time()
  return toc - tic
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
    L,com_lipids,com_chol = trajIO.decompress("comTraj.npz")



    dist1 = edm_two_loop(L[34],com_lipids[34])
    dist2 = edm(L[34],com_lipids[34])

    print(dist1[80][45],dist2[80][45])

    sum = 0.0

    two_loop = 0
    for i in range(5):
        two_loop += time_function(edm_two_loop,L[34],com_lipids[34])
    print('two loop version took %f seconds' % two_loop)

    one_loop = 0
    for i in range(5):
        one_loop += time_function(edm,L[34],com_lipids[34])
    print('one loop version took %f seconds' % one_loop)

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

