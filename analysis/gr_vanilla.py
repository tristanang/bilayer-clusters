import sys
import numpy as np

import trajIO
import constants as c

def gr(com_lipids,com_chol):
	pass

if __name__ == "__main__":
	trajFileName = sys.argv[1]
	Nconf = int(sys.argv[2])
	nlog = int(sys.argv[3])
	topology = sys.argv[4]

	Nchol = trajIO.cholConc(topology)
	N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)

	Nlipids = com_lipids.shape[0]

	com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol)