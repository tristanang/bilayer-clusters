import sys
sys.path.append("../")

from bilayer_clusters import trajIO
from analysis import compress

import numpy


if __name__ == "__main__":
    trajFileName = "../data/temp20"
    chol = "../data/20chol.top"

    Nchol = trajIO.cholConc(chol)
    outfile = "../test/comTraj.npz"

    compress.spitnumpy(trajFileName,chol,100)

    N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,3,100)
    L_new,lipids_new,chol_new = trajIO.decompress(outfile)

    #print(numpy.array_equal(L,L_new))
    #print(numpy.array_equal(com_lipids,lipids_new))
    #print(numpy.array_equal(com_chol,chol_new))

    if numpy.array_equal(L,L_new) and numpy.array_equal(com_lipids,lipids_new) and numpy.array_equal(com_chol,chol_new):
        print("compressLoad test passed")
    else:
        print("compressLoad test failed")

