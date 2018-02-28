import sys
import numpy as np
import pickle

from bilayer_clusters import jenks_clusters
from bilayer_clusters import displacement
from bilayer_clusters import trajIO

if __name__ == "__main__":
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("Needs 3 or 4 arguments")
        print(len(sys.argv))
        exit()

    trajFileName = sys.argv[1]
    Nconf = int(sys.argv[2])
    nlog = int(sys.argv[3])
    topology = sys.argv[-1] 

    if trajIO.rawOrCOM(trajFileName):
        Nchol = trajIO.cholConc(topology)
        N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)

    else:
        L,com_lipids,com_chol = trajIO.decompress(trajFileName)
        Nchol = com_chol.shape[1]

    Nlipids = com_lipids.shape[1]
    Nblock = Nconf//nlog
    
    com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol)

    #parameters
    cluster_sizes = [3,4,5,6,7,8]
    times = list(range(29,38))

    #calculating displacement

    com_lipids = displacement.block_displacement(L,com_lipids)
    com_chol = displacement.block_displacement(L,com_chol)

    #cluster dict
    clusters = {} #cluster[nblock][time][type][layer][cluster_size]
    for block in range(Nblock):
        clusters[block] = {}
        for t in times:
            clusters[block][t] = {}
            clusters[block][t]['lipids'] = {}
            clusters[block][t]['chol'] = {}
            for layer in ['upper','lower']:
                clusters[block][t]['lipids'][layer] = {}
                clusters[block][t]['chol'][layer] = {}

    #running

    for block in range(Nblock):
        for time in times:
            for size in cluster_sizes:
                t = block*nlog + time

                upper_lipids, lower_lipids = trajIO.layering(com_lipids[t])
                upper_chol, lower_chol = trajIO.layering(com_chol[t])
                
                clusters[block][time]['lipids']['upper'][size] = jenks_clusters.clusters(upper_lipids,size)
                clusters[block][time]['lipids']['lower'][size] = jenks_clusters.clusters(lower_lipids,size)
                
                clusters[block][time]['chol']['upper'][size] = jenks_clusters.clusters(upper_chol,size)
                clusters[block][time]['chol']['lower'][size] = jenks_clusters.clusters(lower_chol,size)
            

    output = "clusters.dict"

    pickle.dump(clusters, open(output, "wb" ) )







