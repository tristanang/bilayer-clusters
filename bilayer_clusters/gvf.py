import numpy as np

def gvf(cluster):
    cluster = list(map(lambda x: x[:,3],cluster))

    return gvf_helper(cluster)

def gvf_helper(cluster): #can parallelize for speed
    Nlipids = sum(list(map(lambda x: x.shape[0],cluster)))
    
    pre_mean = sum(list(map(np.sum,cluster)))
    global_mean = pre_mean/Nlipids
    
    class_means = list(map(np.mean,cluster))
    
    SDAM = list(map(lambda x : (x-global_mean)**2, cluster))
    SDAM = sum(list(map(sum, SDAM)))

    for i in range(len(cluster)):
        cluster[i] = (cluster[i] - class_means[i])**2

    SDCM = sum(list(map(np.sum, cluster))) 

    return (SDAM-SDCM)/SDAM
"""
def goodness_of_variance_fit(array, classes):
    # get the break points
    breaks = jenks(array, classes)

    # do the actual classification
    classified = np.array([classify(i, classes) for i in array])

    # max value of zones
    maxz = max(classified)

    # nested list of zone indices
    zone_indices = [[idx for idx, val in enumerate(classified) if zone + 1 == val] for zone in range(maxz)]

    # sum of squared deviations from array mean
    sdam = np.sum((array - array.mean()) ** 2)

    # sorted polygon stats
    array_sort = [np.array([array[index] for index in zone]) for zone in zone_indices]

    # sum of squared deviations of class means
    sdcm = sum([np.sum((classified - classified.mean()) ** 2) for classified in array_sort])

    # goodness of variance fit
    gvf = (sdam - sdcm) / sdam

    return gvf

def classify(value, breaks):
    for i in range(1, len(breaks)):
        if value < breaks[i]:
            return i
    return len(breaks) - 1
"""