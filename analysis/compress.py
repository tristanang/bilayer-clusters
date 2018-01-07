import trajIO
import sys
import constants as c

def spitnumpy(trajFileName, topology, Nconf):
    outfile = "comTraj"
    Nchol = trajIO.cholConc(topology)
    N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)
    
    trajIO.np.savez(outfile,L=L,com_lipids=com_lipids,com_chol=com_chol)

if __name__ == "__main__":
    trajFileName = sys.argv[1]
    Nconf = int(sys.argv[2])
    topology = sys.argv[3]

    spitnumpy(trajFileName, topology, Nconf)

    #print json to file
