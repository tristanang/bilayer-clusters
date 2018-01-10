from bilayer_clusters import jenks
import numpy as np

def clusters_base(arr,k): #array is pre-sorted
    breaks = jenks.jenks(arr,k)
    clusters = [0 for i in range(k)]

    i_start = 0
    i = 0

    for index in breaks[1:-1]:
        i_end = index
        clusters[i] = arr[i_start:i_end+1]
        i += 1
        i_start = i_end+1

    clusters[i] = arr[i_start:]

    return clusters

def clusters(lipids,k,t):
    return

def clusterForm(group,interval,endblock):
    clusters = []
    i_start = 0
    for i in range(1,len(interval)-1):
        i_end = lib.index(list(map(lambda x : x.getS(endblock),group)),interval[i])
        cluster = group[i_start:i_end+1]
        i_start = i_end + 1
        clusters.append(cluster)
    clusters.append(group[i_start:])
    
    return clusters