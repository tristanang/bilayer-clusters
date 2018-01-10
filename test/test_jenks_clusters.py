import json
from bilayer_clusters import jenks_clusters

def test_json():
    data = json.load(open('jenks.json'))
    data.sort()

    cluster1 = data[:442+1]
    cluster2 = data[443:858+1]
    cluster3 = data[859:1246+1]
    cluster4 = data[1247:1604+1]
    cluster5 = data[1605:]

    cluster = [cluster1,cluster2,cluster3,cluster4,cluster5]
    
    clusterlen = list(map(len,cluster))
    assert len(data) == sum(clusterlen)
    #this proves our manual cluster is correct.

    clusterfunc = jenks_clusters.clusters_base(data,5)
    clusterfunclen = list(map(len,clusterfunc))

    assert clusterfunclen == clusterlen
    
    assert clusterfunc == cluster

    return

if __name__ == '__main__':
    test_json()
    print("passed")

#[0, 442, 858, 1246, 1604, 1999]