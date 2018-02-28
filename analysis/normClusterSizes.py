import sys
import pickle
import numpy as np

from bilayer_clusters import trajIO, iter_cluster

if __name__ == '__main__':
    
    trajFileName = sys.argv[1]
    Nconf = 4600
    nlog = 46
    Nblock = Nconf//nlog
    cutoff = 1.15 #anything above 20chol #maybe 1.3 for everything below?

    if trajIO.rawOrCOM(trajFileName):
        Nchol = trajIO.cholConc(topology)
        N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)

        Nlipids = com_lipids.shape[0]
    else:
        L,com_lipids,com_chol = trajIO.decompress(trajFileName)
    
    if True:
        pickle_off = open("clusters.dict","rb")
        clusters = pickle.load(pickle_off)
        #cluster[nblock][time][layer][type][cluster_size]

    #parameters
    cluster_sizes = [3,4,5,6,7,8]
    times = list(range(29,38))

    #initialize output dict
    normSizes = {}
    for block in range(Nblock):
        normSizes[block] = {}
        for t in times:
            normSizes[block][t] = {}
            for particle in ['lipids','chol']:
                normSizes[block][t][particle] = {}
                for layer in ['upper','lower']:
                    normSizes[block][t][particle][layer] = {}
                    for size in cluster_sizes:
                        normSizes[block][t][particle][layer][size] = [0 for i in range(size)]

    weightedNormSizes = normSizes.copy()

    for block in range(Nblock):
        for time in times:
            t = block*nlog + time

            for layer in ['upper','lower']:
                for size in cluster_sizes:
                    for i in range(size):
                        group = len(clusters[block][time]['lipids'][layer][size][i])
                        c = iter_cluster.Cluster(clusters[block][time]['lipids'][layer][size][i],L[t],cutoff)
                        
                        normSizes[block][time]['lipids'][layer][size][i],weightedNormSizes[block][time]['lipids'][layer][size][i] = c.normalize(com_lipids[t])


    output = "normClust.dict"
    output2 = "weightedClust.dict"
    pickle.dump(normSizes, open(output, "wb" ) )
    pickle.dump(weightedNormSizes, open(output, "wb" ) )