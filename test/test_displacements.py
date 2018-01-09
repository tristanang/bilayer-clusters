import sys

from bilayer_clusters import trajIO,displacement
import numpy as np 

def vs02loop():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    noloop = displacement.block_displacement(L,com_lipids)
    twoloop = displacement.block_displacement_loop(L,com_lipids)

    sum = 0.0

    for t in range(52):
        difference = np.linalg.norm(noloop[t] - twoloop[t], ord='fro')
        sum += difference

    if sum < 0.00000001:
        print(str(sum)+" 2loop vs 0 loop passed")
        return True
    else:
        print(str(sum)+" 2loop vs 0 loop failed")
        return False

if __name__ == "__main__":
    vs02loop()