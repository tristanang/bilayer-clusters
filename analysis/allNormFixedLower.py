import sys
import pickle
import numpy as np
from copy import deepcopy

from bilayer_clusters import trajIO
from bilayer_clusters import displacement
from bilayer_clusters import constants as c
from bilayer_clusters import percentages
from bilayer_clusters import density_cluster as dc

def printer(normSizes,weightedSizes,logNorm,logWeighted,Nconf,times,nlog=c.nlog):

    cluster_sizes = c.cluster_sizes
    Nblock = Nconf//nlog

    norm = {}

    for time in c.alltimes:
        norm[time] = {}
        
        for layer in ['upper','lower','together']:
            norm[time][layer] = {}
            
            for size in cluster_sizes:
                norm[time][layer][size] = [[] for i in range(size)]

    weighted = deepcopy(norm)

    for block in range(Nblock):
        times1 = list(normSizes[block].keys())

        for time in times1:
            for layer in ['upper','lower']:
                for size in cluster_sizes:
                    for i in range(size):
                        norm[time][layer][size][i].append(normSizes[block][time][layer][size][i])
                        weighted[time][layer][size][i].append(weightedSizes[block][time][layer][size][i])

    for block in range(Nblock):
        for time in c.alltimes:
            for size in cluster_sizes:
                for i in range(size):
                    norm[time]['together'][size][i] = norm[time]['upper'][size][i] + norm[time]['lower'][size][i]
                    weighted[time]['together'][size][i] = weighted[time]['upper'][size][i] + weighted[time]['lower'][size][i]

    printing = {}
    for kind in ['norm','weighted']:
        printing[kind] = {}
        for size in cluster_sizes:
            printing[kind][size] = [[] for i in range(size)]

    for time in c.alltimes:
        for size in cluster_sizes:
            for i in range(size):
                if norm[time]['together'][size][i]:
                    printing['norm'][size][i].append((np.mean(norm[time]['together'][size][i]))/logNorm[size][i])
                    printing['weighted'][size][i].append((np.mean(weighted[time]['together'][size][i]))/logWeighted[size][i])

    fmt = '%f %f\n'

    for kind in ['norm','weighted']:
        for size in cluster_sizes:
            for i in range(size):
                filename = 'all'+kind+'Size'+'_'+str(size)+'_'+str(i)+".dat"
                f = open(filename,'w')
                for time,data in zip(c.realtimes[:len(printing[kind][size][i])],printing[kind][size][i]):
                    f.write(fmt %(time,data))

                f.close()

if __name__ == '__main__':
    
    trajFileName = sys.argv[1]
    Nconf = int(sys.argv[2])
    nlog = int(sys.argv[3])
    Nblock = Nconf//nlog
    cutoff = 1.3 #anything above 20chol #maybe 1.3 for everything below?

    if trajIO.rawOrCOM(trajFileName):
        Nchol = trajIO.cholConc(topology)
        N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)
        com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol) 

        Nlipids = com_lipids.shape[1]
    else:
        L,com_lipids,com_chol = trajIO.decompress(trajFileName)
        com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol) 
        Nchol = com_lipids.shape[1]

    #parameters
    cluster_sizes = [4]
    times = list(range(1,46))

    if Nchol: com_lipids = np.concatenate((com_lipids,com_chol),axis=1)
    com_lipids = displacement.block_displacement(L,com_lipids)

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

    for block in range(Nblock):
        t = nlog
        time = block*nlog + t

        while time<Nconf:

            normSizes[block][t] = {}

            for layer in ['upper','lower']:
                normSizes[block][t][layer] = {}
                for size in cluster_sizes:
                    normSizes[block][t][layer][size] = [0 for i in range(size)]

            t += nlog
            time += nlog

    weightedNormSizes = deepcopy(normSizes)

    logNorm = {}
    logWeighted = {}
    linearNorm = {}
    linearWeighted = {}

    for size in cluster_sizes:
        logNorm[size] = np.zeros(size)
        logWeighted[size] = np.zeros(size)
        #linearNorm[size] = np.zeros(size)
        #linearWeighted[size] = np.zeros(size)

    #block
    for block in range(Nblock):
        start = block*nlog
        for time in times:
            t = start + time
            #print(t) #progress tracker
            upper,lower = trajIO.layering(com_lipids[t])
            original = {}
            original['upper'] = upper
            original['lower'] = lower

            for layer in ['upper','lower']:
                #clustering
                clusters = percentages.cluster(original[layer],c.percentages['all']['lower'][4])

                for size in cluster_sizes:
                    for i in range(size):
                        Nparticles = len(clusters[i])
                        normSizes[block][time][layer][size][i],weightedNormSizes[block][time][layer][size][i] = dc.mean_cluster_size(clusters[i],L[t],cutoff)
                        
                        if time == 1:
                            alpha,beta = dc.meanRandom(original[layer],L[t],cutoff,Nparticles)
                            logNorm[size][i] += alpha
                            logWeighted[size][i] += beta

    #linear
    #linearNorm = linearWeighted = 

    for block in range(Nblock):
        start = block*nlog
        linear_t = displacement.linear_gen(start,Nconf)

        com_lipids = displacement.linear_displacement(L,com_lipids,start,Nconf)
        
        for time in linear_t:
            print(start,time)
            t = start + time

            upper,lower = trajIO.layering(com_lipids[t])
            original = {}
            original['upper'] = upper
            original['lower'] = lower

            for layer in ['upper','lower']:
                #clustering
                clusters = percentages.cluster(original[layer],c.percentages['all']['lower'][4])

                for size in cluster_sizes:
                    for i in range(size):
                        Nparticles = len(clusters[i])
                        normSizes[block][time][layer][size][i],weightedNormSizes[block][time][layer][size][i] = dc.mean_cluster_size(clusters[i],L[t],cutoff)
                        
    for size in cluster_sizes:
        logNorm[size] /= (Nblock*2)
        logWeighted[size] /= (Nblock*2)
        print(logNorm[size])
        print(logWeighted[size])
        #linearNorm[size] /= (Nblock*2)
        #linearWeighted[size] /= (Nblock*2)

    printer(normSizes,weightedNormSizes,logNorm,logWeighted,Nconf,times,nlog)



    output = "normClust.dict"
    output2 = "weightedClust.dict"

    f = open(output, "wb" )
    pickle.dump(normSizes, f)
    f.close()

    f = open(output2, "wb" )
    pickle.dump(weightedNormSizes, f)
    f.close()
