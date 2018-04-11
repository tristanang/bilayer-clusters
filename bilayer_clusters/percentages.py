import numpy as np 
from bilayer_clusters import constants as c

def split(Nlipids,arr):
    size = len(arr)
    output = [0 for i in range(size)]
    for i in range(size):
        output[i] = round(Nlipids*arr[i])

    difference = Nlipids - sum(output)

    if difference > 0:
        output[-1] += difference

    else:
        output[0] += difference

    assert sum(output) == Nlipids

    return output

def cluster(lipids,percentage):
    Nlipids = len(lipids)
    splits = split(Nlipids,percentage)

    lipids = lipids[lipids[:,3].argsort()]

    output = [0 for i in range(len(splits))]

    output[0] = lipids[:splits[0]]
    end = splits[0]
    for i in range(1,len(splits)):
        output[i] = lipids[end:(end+splits[i])]
        end += splits[i]

    return output



