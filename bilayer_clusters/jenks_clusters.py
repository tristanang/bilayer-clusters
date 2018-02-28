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

def clusters(lipids,k): #removed t dependence
    lipids = lipids[lipids[:,3].argsort()] #weird sorting thing from: https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column

    breaks = jenks.jenks(lipids[:,3],k)
    clusters = [0 for i in range(k)]

    i_start = 0
    i = 0

    for index in breaks[1:-1]:
        i_end = index
        clusters[i] = np.asarray(lipids[i_start:i_end+1])
        i += 1
        i_start = i_end+1

    clusters[i] = np.asarray(lipids[i_start:])

    return clusters
