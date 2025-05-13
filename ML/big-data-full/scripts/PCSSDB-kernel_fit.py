#!/usr/bin/python

#
# This script fits a per-cluster kernel for sum(-ICOHP)
# 
# usage:
# 
# python PCDB-kernel_fit.py <CHEM_COMPOSITION> <CENTRAL_SPECIE> <KERNEL_DIM> <ZETA>
#

# Libraries
import sys
import numpy as np
import pickle as pkl
import quippy

#
# The script
#

# Constants and arguments
BD_DIR = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full"
CHEM_COMPOSITION = sys.argv[1]
DB_DIR = BD_DIR+"/"+CHEM_COMPOSITION+"-PCDB"
CENTRAL_SPECIE = sys.argv[2]
KERNEL_DIM = int(sys.argv[3]) # defines the size of the training set
# https://stackoverflow.com/questions/5018552/how-to-raise-a-numpy-array-to-a-power-corresponding-to-repeated-matrix-multipl/5018638
ZETA = int(sys.argv[4]) # the hyperparametr

# Loading the clusters and their SOAPs
clusters = quippy.AtomsList(DB_DIR+"/"+CENTRAL_SPECIE+".xyz")
ftl = open(DB_DIR+"/"+CENTRAL_SPECIE+"-SOAPS.vec","rb")
soaps = pkl.load(ftl)
ftl.close()
#for i in range(0,10):
#    print clusters[i].params["ID_RUN"], clusters[i].params["SUB_STEP"], clusters[i].params["ID_CENTRAL_ATOM"], np.dot(soaps[0], soaps[i])

# Building the kernel matrix
KM = np.zeros((KERNEL_DIM,KERNEL_DIM))
for i in range(KERNEL_DIM):
    for j in range(KERNEL_DIM):
        KM[i][j] = np.dot(soaps[i], soaps[j])

# Now the variance of the sum(-ICOHP) values
ICOHPs_TRAINING = []
for i in range(KERNEL_DIM):
    ICOHPs_TRAINING.append(float(clusters[i].params["SUM_ICOHP"]))
VAR = np.var(ICOHPs_TRAINING)

# Now setting up equation (3) in the supplementary material of  [Nature Communications 9, 4501 (2018)]
# Firstly the inverse matrix ... 
IVAR = VAR * np.identity(KERNEL_DIM)
KMM = np.matrix(KM)
IKM = np.linalg.inv((KMM**ZETA) + np.matrix(IVAR))
#IKM = np.linalg.inv(np.linalg.matrix_power(KM,ZETA) + IVAR)
# ... and then weights (alpha_i)
alpha = [0.0] * KERNEL_DIM
for i in range(KERNEL_DIM):
    for j in range(KERNEL_DIM):
        alpha[i] += IKM.item((i,j)) * ICOHPs_TRAINING[j]
        #alpha[i] += IKM[i][j] * ICOHPs_TRAINING[j]

# And finally the predicted sum(-ICOHP) values for the test set
ICOHPs_TEST = [0.0] * (len(clusters) - KERNEL_DIM)
SUMDIFF = 0.0
for it in range(KERNEL_DIM, len(clusters)):
    for i in range(KERNEL_DIM):
        ICOHPs_TEST[it-KERNEL_DIM] += alpha[i]*(np.dot(soaps[it], soaps[i])**ZETA)
    # checking the root-mean square error
    print it-KERNEL_DIM, it, ICOHPs_TEST[it-KERNEL_DIM], clusters[it].params["SUM_ICOHP"]
    SUMDIFF += (ICOHPs_TEST[it-KERNEL_DIM] - clusters[it].params["SUM_ICOHP"])**2
RMSD = np.sqrt(SUMDIFF/len(ICOHPs_TEST))
print RMSD

