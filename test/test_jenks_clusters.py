import json
from bilayer_clusters import jenks_clusters, trajIO, displacement

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

def test_sort():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    com_lipids = displacement.block_displacement(L,com_lipids)
    Nlipids = com_lipids.shape[1]
    Nconf = com_lipids.shape[0]

    for t in range(Nconf):
        lipids = com_lipids[t][com_lipids[t][:,3].argsort()]
        curr = 0
        for i in range(Nlipids):
            if curr > lipids[i,3]:
                print("sorting messed up")
                raise ValueError
            else:
                curr = lipids[i,3]

    return True

def test_clusters(): #not a complete testfunction but the function should nonetheless be accurate
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    com_lipids = displacement.block_displacement(L,com_lipids)
    Nlipids = com_lipids.shape[1]

    clusters = jenks_clusters.clusters(com_lipids[40],4)

    return clusters

if __name__ == '__main__':
    test_json()
    test_sort()
    test_clusters()
    print("passed")

#[0, 442, 858, 1246, 1604, 1999]