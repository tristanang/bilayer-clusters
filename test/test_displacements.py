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

def s_test():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    control = displacement.block_displacement(L,com_lipids)
    control = control[:,:,3]

    test = com_lipids

    test = displacement.s(L,com_lipids,0,list(range(1,46)))
    test = displacement.s(L,test,46,np.asarray(list(range(1,46)))+46)
    test = displacement.s(L,test,92,[93,94,95,96,97,98,99])
    test = test[:,:,3]

    for t in range(100):
        boo = (test[t] == control[t])
        print(boo.all())

def linear_test():
    file = "comTraj.npz"
    L,com_lipids,com_chol = trajIO.decompress(file)

    control = displacement.block_displacement(L,com_lipids)
    control1 = control[:,:,3]

    test = displacement.linear_displacement(L,control,0,100)
    test1 = test[:,:,3]
    print(test1.shape)

    for t in range(46):
        boo = (test1[t] == control1[t])
        print(boo.all())

    print("hi")

    for t in [46,92]:
        boo = (test1[t] != control1[t])
        print(boo.all())


if __name__ == "__main__":
    #vs1_2loop()
    #vs0_2loop()
    #speed benchmark

    #timings()
    #assertion = displacement.linear_gen(4554,4600)
    
    #s_test()
    linear_test()

