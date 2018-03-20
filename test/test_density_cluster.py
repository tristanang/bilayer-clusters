from bilayer_clusters import density_cluster as dc
import sys

from bilayer_clusters import trajIO

def time_function(f, *args):
  """
  Call a function f with args and return the time (in seconds) that it took to execute.
  """
  import time
  tic = time.time()
  f(*args)
  toc = time.time()
  return toc - tic

if __name__ == '__main__':
    L, com_lipids, com_chol = trajIO.decompress("comTraj.npz")
    import pickle

    pickle_off = open("clusters.dict","rb")
    clusters = pickle.load(pickle_off)

    X = clusters[0][28]['lipid'][0][0]

    db_time = 0
    for i in range(10):
        db_time += time_function(dc.dbscan_wrapper, X,L[0],1.15)
    print('DBSCAN version took %f seconds' % db_time)

    hierarchy_time = 0
    for i in range(10):
        hierarchy_time += time_function(dc.hierarchy_wrapper, X,L[0],1.15)
    print('Hierarchy version took %f seconds' % hierarchy_time)

