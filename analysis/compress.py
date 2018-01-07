import trajIO
import sys
import json
import constants as c

def spitjson(trajFileName, topology, Nconf):
	Nchol = trajIO.cholConc
	N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)
	#json dump

if __name__ == "__main__":
	trajFileName = sys.argv[1]
	Nconf = int(sys.argv[2])
	nlog = int(sys.argv[3])
	topology = sys.argv[4]

	#print json to file
