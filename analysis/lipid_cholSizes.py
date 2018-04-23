from bilayer_clusters import trajIO
from bilayer_clusters import jenks_clusters
from bilayer_clusters import displacement

from copy import deepcopy
import numpy as np
import pickle

def sizet(thing,size,arr):
    colors = ['b', 'c', 'y', 'm', 'r','k','b','c']
    times = [29,30,31,32,33,34,35,36,37]
    
    lst []

    for i in range(size):
        y = []
    
        for time in times:
            y.append(np.mean(arr[thing][time][size][i]))
        
        y = np.asarray(y)
        lst.append(np.mean(y))

    return lst
    
if __name__ == '__main__':
    import sys
    sys.stdout = open('sizes.txt','wt')
        
    Nblock = 100

    pickle_off1 = open("clusters.dict","rb")
    pickle_off2 = open("clusters15.dict","rb")

    cluster_sizes = [2,3,4,5]
    clusters = pickle.load(pickle_off1)
    clusters15 = pickle.load(pickle_off2)

    times = [29,30,31,32,33,34,35,36,37]

    c_size = {}
    c_size['lipids'] = {}

    for thing in ['lipids']:
        for time in times:
            c_size[thing][time] = {}
            for size in [2,3,4,5]:
                c_size[thing][time][size] = [[] for i in range(size)]

    for block in range(Nblock):
        times1 = list(clusters[block].keys())

        for time in times:
            for size in [3,4]:
                for layer in ['upper','lower']:
                    temp = list(map(lambda x: len(x),clusters[block][time]['lipids'][layer][size]))
                    temp = np.asarray(temp)
                    summ = np.sum(temp)
                    temp = temp / summ
                    
                    for i in range(size):
                        c_size['lipids'][time][size][i].append(temp[i])


    for block in range(Nblock):
        times1 = list(clusters[block].keys())
        
        for time in times:
            for size in [2,5]:
                for layer in ['upper','lower']:
                    temp = list(map(lambda x: len(x),clusters15[block][time]['lipids'][layer][size]))
                    temp = np.asarray(temp)
                    summ = np.sum(temp)
                    temp = temp / summ
                    
                    for i in range(size):
                        c_size['lipids'][time][size][i].append(temp[i])

    for size in cluster_sizes:
        sizet('lipids',size,c_size)
        #sizet('chol',size)
        #plt.show()
        print("----------"+str(size)+"-----------")