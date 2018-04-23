from bilayer_clusters import trajIO
from bilayer_clusters import constants as c
import sys

import numpy as np


def spitnumpy(trajFileName1, trajFileName2, Nconf):
    outfile = "tailTraj"
    NDIM = 3

    L = np.zeros((Nconf,NDIM), dtype=np.double)
    com_lipids = np.zeros((Nconf,Nlipids,NDIM+1), dtype=np.double)
    com_chol = np.zeros((Nconf,Nchol,NDIM+1), dtype=np.double)

    trajFile1 = open(trajFileName1,'r')
    Nlipids = int(trajFile1.readline().strip())
    print(Nlipids)

    for t in range(Nconf):
        time = trajFile1.readline()
        assert t == time

        line = trajFile1.readline().split()
        for k in range(NDIM):
            L[t,k] = float(line[k])

        for i in range(Nlipids):
            line = trajFile1.readline().split()
            for k in range(NDIM+1):
                com_lipids[t][i][k] = line[k]

    trajFile1.close()

    trajFile2 = open(trajFileName2, 'r')
    Nchol = int(trajFile2.readline().strip())
    print(Nchol)

    for t in range(Nconf):
        time = trajFile.readline()
        assert t == time

        for i in range(Nchol):
            line = trajFile.readline().split()
            for k in range(NDIM+1):
                com_lipids[t][i][k] = line[k]


    com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol)    

    trajIO.np.savez(outfile,L=L,com_lipids=com_lipids,com_chol=com_chol)

if __name__ == "__main__":
    trajFileName1 = sys.argv[1]
    trajFileName2 = sys.argv[2]
    Nconf = int(sys.argv[3])

    spitnumpy(trajFileName1, trajFileName2, Nconf)
