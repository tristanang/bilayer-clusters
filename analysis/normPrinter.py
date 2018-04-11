import pickle
import numpy as np
from copy import deepcopy
import sys   

from bilayer_clusters import constants as c

pickle_off1 = open("normClust.dict","rb")
pickle_off2 = open("weightedClust.dict","rb")

normSizes = pickle.load(pickle_off1)
weightedSizes = pickle.load(pickle_off2)

cluster_sizes = c.cluster_sizes
times = list(range(1,46))

Nconf = int(sys.argv[1])
nlog = int(sys.argv[2])
Nblock = Nconf//nlog

norm = {}

for time in times:
    norm[time] = {}
    
    for layer in ['upper','lower','together']:
        norm[time][layer] = {}
        
        for size in cluster_sizes:
            norm[time][layer][size] = [[] for i in range(size)]
            
weighted = deepcopy(norm)

for block in range(Nblock):
    for time in times:
        for layer in ['upper','lower']:
            for size in cluster_sizes:
                for i in range(size):
                    norm[time][layer][size][i].append(normSizes[block][time][layer][size][i])
                    weighted[time][layer][size][i].append(weightedSizes[block][time][layer][size][i])

for block in range(Nblock):
    for time in times:
        for size in cluster_sizes:
            for i in range(size):
                norm[time]['together'][size][i] = norm[time]['upper'][size][i] + norm[time]['lower'][size][i]
                weighted[time]['together'][size][i] = weighted[time]['upper'][size][i] + weighted[time]['lower'][size][i]

printing = {}
for kind in ['norm','weighted']:
    printing[kind] = {}
    for size in cluster_sizes:
        printing[kind][size] = [[] for i in range(size)]

for time in times:
    for size in cluster_sizes:
        for i in range(size):
            printing['norm'][size][i].append(np.mean(norm[time]['together'][size][i]))
            printing['weighted'][size][i].append(np.mean(weighted[time]['together'][size][i]))

fmt = '%d %f\n'

for kind in ['norm','weighted']:
    for size in cluster_sizes:
        for i in range(size):
            filename = kind+'Size'+'_'+str(size)+'_'+str(i)+".dat"
            f = open(filename,'w')
            for time,data in zip(times,printing[kind][size][i]):
                f.write(fmt %(time,data))

            f.close()
