import sys

from bilayer_clusters import trajIO,displacement
import numpy as np 

def time_function(f, *args):
    """
    Call a function f with args and return the time (in seconds) that it took to execute.
    """
    import time
    tic = time.time()
    f(*args)
    toc = time.time()
    return toc - tic

def vs1_2loop():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    noloop = displacement.block_displacement_one(L,com_lipids)
    twoloop = displacement.block_displacement_loop(L,com_lipids)

    noloop = noloop[:,:,3]
    twoloop = twoloop[:,:,3]

    for t in range(100):
        boo = (noloop[t] == twoloop[t])
        print(boo.all())
                


def vs0_2loop():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    noloop = displacement.block_displacement_no(L,com_lipids)
    twoloop = displacement.block_displacement_loop(L,com_lipids)

    noloop = noloop[:,:,3]
    twoloop = twoloop[:,:,3]

    for t in range(100):
        boo = (noloop[t] == twoloop[t])
        print(boo.all())
    """
    for t in range(100):
        difference = np.linalg.norm(noloop[t] - twoloop[t], ord='fro')
        sum += difference

    if sum < 0.00000001:
        print(str(sum)+" 2loop vs 0 loop passed")
        return True
    else:
        print(str(sum)+" 2loop vs 0 loop failed")
        return False
    """
def timings():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)
    two_loop_time = one_loop_time = no_loop_time = 0

    for t in range(10):
        #two_loop_time += time_function(displacement.block_displacement_loop, L,com_lipids)

        one_loop_time += time_function(displacement.block_displacement_one, L,com_lipids)
    
        no_loop_time += time_function(displacement.block_displacement_no, L,com_lipids)
    
    print('Two loop version took %f seconds' % two_loop_time)
    print('One loop version took %f seconds' % one_loop_time)
    print('No loop version took %f seconds' % no_loop_time)

if __name__ == "__main__":
    #vs1_2loop()
    #vs0_2loop()
    #speed benchmark

    timings()