from jenks import jenks
from bisect import bisect_left
import numpy as np

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def binary_search(arr,x,lo=0,hi=None):
    if not hi:
        hi = len(arr)

    found = False

    while lo<=hi and not found:
        mid = (lo + hi)//2
        if arr[mid] == x:
            return mid
        elif arr[mid]>x:
            hi = mid-1
        elif arr[mid]<x:
            lo = mid+1
        else:
            raise ValueError

    return ValueError

def clusters_base(arr,k): #array is pre-sorted
    breaks = jenks(arr,k)
    print(breaks)
    clusters = []

    breaks = list(map(lambda x : binary_search(breaks,x),breaks))

    return(breaks)

import numpy as np
arr = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,100.0]
arr.sort()
print(arr)
print(clusters_base(arr,3))


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