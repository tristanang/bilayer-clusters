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

    upper = 0
    lower = 1
    
    com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol)

    #parameters
    Ncluster = 4
    #times = [10,19,28,37,46]
    times = [28,37,46]

    #calculating displacement

    com_lipids = displacement.block_displacement(L,com_lipids)

    #cluster forming
    clusters = {} #cluster[t][time]['lipid'][upper][cluster number]
    for t in range(Nblock):
        clusters[t] = {}
        for time in times:
            clusters[t][time] = {}
            clusters[t][time]['lipid'] = [0,1] 
            clusters[t][time]['chol'] = [0,1]

    for start in range(Nblock):
        for time in times:
            t = start + time

            upper_lipids, lower_lipids = trajIO.layering(com_lipids[t])
            upper_chol, lower_chol = trajIO.layering(com_chol[t])

            cluster_u = jenks_clusters.clusters(upper_lipids,Ncluster) 
            cluster_l = jenks_clusters.clusters(lower_lipids,Ncluster)

            for n in range(Ncluster):
                cluster_u[n] = cluster_u[n][:,:2]
                cluster_l[n] = cluster_l[n][:,:2]

            clusters[start][time]['lipid'][upper] = cluster_u
            clusters[start][time]['lipid'][lower] = cluster_l
            clusters[start][time]['chol'][upper] = upper_chol[:,:2]
            clusters[start][time]['chol'][lower] = lower_chol[:,:2]

    output = "clusters.dict"

    pickle.dump(clusters, open(output, "wb" ) )







