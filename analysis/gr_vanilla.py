import sys
import numpy as np

from bilayer_clusters import trajIO
from bilayer_clusters import constants as c
from bilayer_clusters import euclideanDist

def gr(L,shape,matrix1,matrix2=None,edm=euclideanDist.edm):
    flag = False

    if not np.any(matrix2):
        matrix2 = matrix1
        flag = True

    gr = np.zeros(shape,int)

    dists = edm(L,matrix1,matrix2)
    dists = dists//c.DR
    dists = dists.astype(int)

    unique, counts = np.unique(dists, return_counts=True)
    hsh = dict(zip(unique,counts))

    for k,v in hsh.items():
        gr[k] += v

    if flag:
        gr[0] -= matrix1.shape[0]

    return gr

def norm(L,lipid_pair,cross_pair=None,chol_pair=None):
    L_ave = [np.mean(L[:,0]),np.mean(L[:,1])]
    
    lipid_norm = lipid_pair * np.pi * c.DR / 2.
    lipid_norm = L_ave[0]*L_ave[1]/(4*lipid_norm)
    
    cross_norm = cross_pair * np.pi * c.DR /2.
    cross_norm = L_ave[0]*L_ave[1]/(4*cross_norm)

    chol_norm = chol_pair * np.pi * c.DR /2.
    chol_norm = L_ave[0]*L_ave[1]/(4*chol_norm)

    return lipid_norm,cross_norm,chol_norm

if __name__ == "__main__":
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("Needs 3 or 4 arguments")
        print(len(sys.argv))
        exit()

    trajFileName = sys.argv[1]
    Nconf = int(sys.argv[2])
    nlog = int(sys.argv[3])
    topology = sys.argv[-1] 

    if trajIO.rawOrCOM(trajFileName):
        Nchol = trajIO.cholConc(topology)
        N,L,com_lipids,com_chol = trajIO.processTrajCOM(trajFileName,Nchol,c.NDIM,Nconf)

        Nlipids = com_lipids.shape[0]
    else:
        L,com_lipids,com_chol = trajIO.decompress(trajFileName)

    com_lipids,com_chol = trajIO.translateZ(com_lipids,com_chol)

    ##data arrays

    lipid_gr = np.zeros([int(L[0][0]/c.DR)],int)
    cross_gr =  np.zeros([int(L[0][0]/c.DR)],int)
    chol_gr = np.zeros([int(L[0][0]/c.DR)],int)
    
    lipid_pair = 0
    cross_pair = 0
    chol_pair = 0

    for t in range(Nconf):
        if t%nlog == 0:
            upper_lipids, lower_lipids = trajIO.layering(com_lipids[t])

            upper_chol, lower_chol = trajIO.layering(com_chol[t])

            ulc, llc = upper_lipids.shape[0], lower_lipids.shape[0]

            ucc, lcc = upper_chol.shape[0], lower_chol.shape[0]

            #lipid
            lipid_gr += gr(L[t],lipid_gr.shape,upper_lipids)
            lipid_gr += gr(L[t],lipid_gr.shape,lower_lipids)

            lipid_pair += ulc*ulc - ulc
            lipid_pair += llc*llc - llc

            #cross
            cross_gr += gr(L[t],cross_gr.shape,upper_lipids,upper_chol)
            cross_gr += gr(L[t],cross_gr.shape,lower_lipids,lower_chol)

            cross_pair += ulc*ucc
            cross_pair += llc*lcc

            #chol
            chol_gr += gr(L[t],chol_gr.shape,upper_chol)
            chol_gr += gr(L[t],chol_gr.shape,lower_chol)

            chol_pair += ucc*ucc - ucc
            chol_pair += lcc*lcc - lcc


    lipid_norm,cross_norm,chol_norm = norm(L,lipid_pair,cross_pair,chol_pair)

    lipidFile = open("gr-com-lipids.dat",'w')

    for step in range(len(lipid_gr)//2):
        lipidFile.write(str((step+0.5)*c.DR)+" "+str(lipid_norm*lipid_gr[step]/((step+0.5)*c.DR))+"\n")
    lipidFile.close()

    cholFile = open("gr-com-cross.dat",'w')

    for step in range(len(cross_gr)//2):
        cholFile.write(str((step+0.5)*c.DR)+" "+str(cross_norm*cross_gr[step]/((step+0.5)*c.DR))+"\n")       
    cholFile.close()

    cholesterol = open("gr-com-chol.dat",'w')
    for step in range(len(chol_gr)//2):
        cholesterol.write(str((step+0.5)*c.DR)+" "+str(chol_norm*chol_gr[step]/((step+0.5)*c.DR))+"\n")       
    cholesterol.close()






    