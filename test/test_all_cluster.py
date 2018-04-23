from bilayer_clusters import constants as c
from bilayer_clusters import iter_cluster as iter
from bilayer_clusters import density_cluster as dc 
from bilayer_clusters import trajIO
from bilayer_clusters import displacement
from bilayer_clusters import percentages
from bilayer_clusters import euclideanDist

import numpy as np

def main():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)
    com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol)

    com_lipids = displacement.block_displacement(L,com_lipids)
    com_chol = displacement.block_displacement(L,com_chol)
    
    t=28
    lipids = com_lipids[t]
    chol = com_chol[t]

    lipids,trash = trajIO.layering(lipids)
    chol,trash = trajIO.layering(chol)

    total = np.concatenate((lipids, chol), axis=0)
    total1 = iter.combine(lipids,chol)



    cluster = percentages.cluster(total,[0.25,0.25,0.25,0.25])
    cluster1 = percentages.cluster(total1,[0.25,0.25,0.25,0.25])

    #edm = euclideanDist.edm(L[t],cluster[0])
    #edm1 = euclideanDist.edm(L[t],cluster1[0])
    #print(np.array_equiv(edm,edm1))


    cutoff = 1.15

    labels1 = dc.dbscan_wrapper(cluster[0],L[t],cutoff)
    labels2 = iter.cluster_labels('upper',L[t],cluster1[0])

    return labels1,labels2






