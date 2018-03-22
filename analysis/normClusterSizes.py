import sys
import pickle
import numpy as np
from copy import deepcopy

from bilayer_clusters import trajIO
from bilayer_clusters import density_cluster as dc

if __name__ == '__main__':
    
    trajFileName = sys.argv[1]
    Nconf = int(sys.argv[2])
    nlog = int(sys.argv[3])
    Nblock = Nconf//nlog
    cutoff = 1.15 #anything above 20chol #maybe 1.3 for everything below?

    if trajIO.rawOrCOM(trajFileName):
        Nchol = trajIO.cholConc(topology)
        N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)
        

        Nlipids = com_lipids.shape[0]
    else:
        L,com_lipids,com_chol = trajIO.decompress(trajFileName)
        
    com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol) 
    
    if True:
        pickle_off = open("clusters.dict","rb")
        clusters = pickle.load(pickle_off)
        #cluster[nblock][time][layer][type][cluster_size]

    #parameters
    del com_chol
    cluster_sizes = [3,4]
    times = list(range(1,46))

    #initialize output dict
    normSizes = {}
    for block in range(Nblock):
        normSizes[block] = {}
        for t in times:
            normSizes[block][t] = {}
        
            for layer in ['upper','lower']:
                normSizes[block][t][layer] = {}
                for size in cluster_sizes:
                    normSizes[block][t][layer][size] = [0 for i in range(size)]

    weightedNormSizes = deepcopy(normSizes)

    for block in range(Nblock):
        start = block*nlog
        upper, lower = trajIO.layering(com_lipids[start])
        original = {}
        original['upper'] = upper
        original['lower'] = lower
        print(len(upper),len(lower))
        for time in times:
            t = start + time
            print(t) #progress tracker
            
            for layer in ['upper','lower']:
                for size in cluster_sizes:
                    for i in range(size):
                        index = clusters[block][time]['lipids'][layer][size][i]
                        normSizes[block][time][layer][size][i],weightedNormSizes[block][time][layer][size][i] = dc.normSize(original[layer],index,L[t],cutoff)

    output = "normClust.dict"
    output2 = "weightedClust.dict"

    f = open(output, "wb" )
    pickle.dump(normSizes, f)
    f.close()

    f = open(output2, "wb" )
    pickle.dump(weightedNormSizes, f)
    f.close()