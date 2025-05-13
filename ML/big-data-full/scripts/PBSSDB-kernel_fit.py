#!/usr/bin/python

#
# This script fits a per-bond (PB) single-SOAP (SS) kernel for -ICOHP values
# 
# usage:
# 
# python kernel_fit.py <CROSSED> <CHEM_COMPOSITION_TRAIN> [<CHEM_COMPOSITION_TEST>] <SPECIE_A> <SPECIE_B> <KERNEL_DIM> <ZETA> <SIGMA> <REGP> <TEST_SIZE> <DIR_TO_WRITE> <SHOW_PLOT> <SAVE_PLOT>
#

# Libraries
import sys
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import math
from random import shuffle

#
# The script
#

# Constants and arguments
BIGDATA_DIR = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full"
CROSS = sys.argv[1] == "True"
if CROSS: # use one database for testing and another one for training
    iarg = 1
    CHEM_COMPOSITION_TRAIN = sys.argv[2]
    CHEM_COMPOSITION_TEST = sys.argv[3]
else: # training and testing sets from the same database
    iarg = 0
    CHEM_COMPOSITION_TRAIN = sys.argv[2]
    CHEM_COMPOSITION_TEST = sys.argv[2]
SPECIE_A = sys.argv[3+iarg]
SPECIE_B = sys.argv[4+iarg]
KERNEL_DIM = int(sys.argv[5+iarg]) # defines the size of the training set
# https://stackoverflow.com/questions/5018552/how-to-raise-a-numpy-array-to-a-power-corresponding-to-repeated-matrix-multipl/5018638
if sys.argv[6+iarg] == "SE" or sys.argv[6+iarg] == "OU":
    ZETA = sys.argv[6+iarg]
else:
    ZETA = float(sys.argv[6+iarg]) # the hyperparameter
SIGMA = float(sys.argv[7+iarg]) # I think this is the characteristic length-scale of the process
REGP = float(sys.argv[8+iarg]) # A hyperparameter proposed by me
TEST_SIZE = int(sys.argv[9+iarg]) # the testing set size
if len(sys.argv) > 10+iarg:
    DIR_TO_WRITE = sys.argv[10+iarg]
    SHOW_PLOT = sys.argv[11+iarg] == "True" # plot the chart?
    SAVE_PLOT = sys.argv[12+iarg] == "True" # plot the chart?
else:
    DIR_TO_WRITE = ""
    SHOW_PLOT = False
    SAVE_PLOT = False

# Loading the bonds of the training and testing sets and their respective SOAPs
DB_DIR = BIGDATA_DIR+"/"+CHEM_COMPOSITION_TRAIN+"-PBSSDB"
fileobj = open(DB_DIR+"/"+SPECIE_A+"-"+SPECIE_B+".bnd")
bondsff = fileobj.readlines()
fileobj.close()
# shuffling the database before using it
shuffle(bondsff)
if not CROSS and KERNEL_DIM+TEST_SIZE > len(bondsff):
    sys.exit("The sum <KERNEL_DIM>+<TEST_SIZE> is greater than the total number of bonds available in the database!!!")
elif KERNEL_DIM > len(bondsff):
    sys.exit("The <KERNEL_DIM> is greater than the total number of bonds available in the training database!!!")
BONDS_TRAIN = bondsff[:KERNEL_DIM]
ftl = open(DB_DIR+"/"+SPECIE_A+"-"+SPECIE_B+"-SOAPS.vec","rb")
SOAPS_TRAIN = pkl.load(ftl)
ftl.close()
if CROSS:
    DB_DIR = BIGDATA_DIR+"/"+CHEM_COMPOSITION_TEST+"-PBSSDB"
    fileobj = open(DB_DIR+"/"+SPECIE_A+"-"+SPECIE_B+".bnd")
    bondsff = fileobj.readlines()
    fileobj.close()
    # shuffling the database before using it
    shuffle(bondsff)
    if TEST_SIZE > len(bondsff):
        sys.exit("The <TEST_SIZE> is greater than the total number of bonds available in the testing database!!!")
BONDS_TEST = bondsff[len(bondsff)-TEST_SIZE:]
if CROSS:
    ftl = open(DB_DIR+"/"+SPECIE_A+"-"+SPECIE_B+"-SOAPS.vec","rb")
    SOAPS_TEST = pkl.load(ftl)
    ftl.close()
else:
    SOAPS_TEST = SOAPS_TRAIN # reference? or duplicated?
	
# Building the kernel matrix
KM = np.zeros((KERNEL_DIM,KERNEL_DIM))
for i in range(KERNEL_DIM):
    bi_dist = float(BONDS_TRAIN[i].split()[4])
    soapAi = SOAPS_TRAIN[str(BONDS_TRAIN[i].split()[0])+"-"+str(BONDS_TRAIN[i].split()[1])+"-"+str(BONDS_TRAIN[i].split()[2])]
    soapBi = SOAPS_TRAIN[str(BONDS_TRAIN[i].split()[0])+"-"+str(BONDS_TRAIN[i].split()[1])+"-"+str(BONDS_TRAIN[i].split()[3])]
    for j in range(KERNEL_DIM):
        bj_dist = float(BONDS_TRAIN[j].split()[4])
        soapAj = SOAPS_TRAIN[str(BONDS_TRAIN[j].split()[0])+"-"+str(BONDS_TRAIN[j].split()[1])+"-"+str(BONDS_TRAIN[j].split()[2])]
        soapBj = SOAPS_TRAIN[str(BONDS_TRAIN[j].split()[0])+"-"+str(BONDS_TRAIN[j].split()[1])+"-"+str(BONDS_TRAIN[j].split()[3])]
        dpAAij = np.dot(soapAi,soapAj)
        dpBBij = np.dot(soapBi,soapBj)
        dpABij = np.dot(soapAi,soapBj)
        dpBAij = np.dot(soapBi,soapAj)
        if sys.argv[6+iarg] == "SE":
            KM[i][j] = np.exp((-(bi_dist-bj_dist)**2)/(2*SIGMA**2)) # Squared Exponential (SE)
        elif sys.argv[6+iarg] == "OU":
            KM[i][j] = np.exp(-abs(bi_dist-bj_dist)/SIGMA) # Ornstein-Uhlenbeck (OU)
        else:
            KM[i][j] = np.exp((-(bi_dist-bj_dist)**2)/(2*SIGMA**2))*(0.25*(dpAAij+dpBBij+dpABij+dpBAij))**ZETA

# Now the variance and standard deviation of the -ICOHP values in the testing set
ICOHPs_TRAIN = []
for i in range(KERNEL_DIM):
    ICOHPs_TRAIN.append(float(BONDS_TRAIN[i].split()[5]))
STD = np.std(ICOHPs_TRAIN)
VAR = STD**2

# Now setting up equation (3) in the supplementary material of [Nature Communications 9, 4501 (2018)]
# Firstly the inverse matrix ... 
#IVAR = VAR * np.identity(KERNEL_DIM)
# on 09/07/2019 Prof. Gabor said to do this
IREG = STD * REGP * np.identity(KERNEL_DIM)
# Note: using KMM and np.linalg.inv below doesn"t matter, see
# ~/UFSCar/MG-NMR/ML/big-data-full/scripts/Zr45Cu45Al10_tests/PB/3/Zr45Cu45Al10-Zr-Al-1500-2-0.5
# ~/UFSCar/MG-NMR/ML/big-data-full/scripts/Zr45Cu45Al10_tests/PB/3/Zr45Cu45Al10-Zr-Al-1500-2-0.5_TEST_np_matrix_power
# Here, I"m trying to use ZETA according to [Nature Communications 9, 4501 (2018)]
#KMM = np.matrix(KM)
#IKM = np.linalg.inv((KMM**ZETA) + np.matrix(IVAR))
#IKM = np.linalg.inv(np.linalg.matrix_power(KM,ZETA) + IVAR)
# Here I will use ZETA latter, as suggested by Prof. Gabor on 05/07/2019
#IKM = np.linalg.inv(KM + IVAR)
# on 09/07/2019 Prof. Gabor said to do this
IKM = np.linalg.inv(KM + IREG)
# ... and then weights (alpha_i)
alpha = [0.0] * KERNEL_DIM
for i in range(KERNEL_DIM):
    for j in range(KERNEL_DIM):
        alpha[i] += IKM.item((i,j)) * ICOHPs_TRAIN[j]
        #alpha[i] += IKM[i][j] * ICOHPs_TRAIN[j]

# And finally the predicted sum(-ICOHP) values for the test set
if DIR_TO_WRITE != "":
    if CROSS:
        FILES_PREFIX = CHEM_COMPOSITION_TRAIN+"-"+CHEM_COMPOSITION_TEST+"-"+SPECIE_A+"-"+SPECIE_B+"-"+str(KERNEL_DIM)+"-"+str(ZETA)+"-"+str(SIGMA)+"-"+str(REGP)+"-"+str(TEST_SIZE)
    else:
        FILES_PREFIX = CHEM_COMPOSITION_TRAIN+"-"+SPECIE_A+"-"+SPECIE_B+"-"+str(KERNEL_DIM)+"-"+str(ZETA)+"-"+str(SIGMA)+"-"+str(REGP)+"-"+str(TEST_SIZE)
    dattw = open(DIR_TO_WRITE+"/"+FILES_PREFIX+".dat","w")
    ftw = open(DIR_TO_WRITE+"/"+FILES_PREFIX+".info","w")
    dtpx = []
    dtpy = []
    higher = 0.0

# Now testing
ICOHPs_TEST = [0.0] * TEST_SIZE
SUMDIFF = 0.0
# The testing set will be the lattest BONDS_TEST in the array
tit = 0 # the true index in the array ICOHPs_TEST
for it in range(TEST_SIZE):
    bit_dist = float(BONDS_TEST[it].split()[4])
    soapAit = SOAPS_TEST[str(BONDS_TEST[it].split()[0])+"-"+str(BONDS_TEST[it].split()[1])+"-"+str(BONDS_TEST[it].split()[2])]
    soapBit = SOAPS_TEST[str(BONDS_TEST[it].split()[0])+"-"+str(BONDS_TEST[it].split()[1])+"-"+str(BONDS_TEST[it].split()[3])]
    for i in range(KERNEL_DIM):
        bi_dist = float(BONDS_TRAIN[i].split()[4])
        soapAi = SOAPS_TRAIN[str(BONDS_TRAIN[i].split()[0])+"-"+str(BONDS_TRAIN[i].split()[1])+"-"+str(BONDS_TRAIN[i].split()[2])]
        soapBi = SOAPS_TRAIN[str(BONDS_TRAIN[i].split()[0])+"-"+str(BONDS_TRAIN[i].split()[1])+"-"+str(BONDS_TRAIN[i].split()[3])]
        dpAAiti = np.dot(soapAit,soapAi)
        dpBBiti = np.dot(soapBit,soapBi)
        dpABiti = np.dot(soapAit,soapBi)
        dpBAiti = np.dot(soapBit,soapAi)
        # Here I"m using ZETA as suggested by Prof. Gabor on 05/07/2019
        if sys.argv[6+iarg] == "SE":
            KERNEL = np.exp((-(bit_dist-bi_dist)**2)/(2*SIGMA**2)) # Squared Exponential (SE)
        elif sys.argv[6+iarg] == "OU":
            KERNEL = np.exp(-abs(bit_dist-bi_dist)/SIGMA) # Ornstein-Uhlenbeck (OU)
        else:
            KERNEL = np.exp((-(bit_dist-bi_dist)**2)/(2*SIGMA**2))*(0.25*(dpAAiti+dpBBiti+dpABiti+dpBAiti))**ZETA
        # Here, I"m trying to use ZETA according to [Nature Communications 9, 4501 (2018)]
        # ICOHPs_TEST[tit] += alpha[i]*(KERNEL**ZETA)
        # Here I already used ZETA before, as suggested by Prof. Gabor on 05/07/2019
        ICOHPs_TEST[tit] += alpha[i]*KERNEL
    SUMDIFF += (ICOHPs_TEST[tit] - float(BONDS_TEST[it].split()[5]))**2
    if DIR_TO_WRITE != "":
        dattw.write("%10f %10f\n" % (ICOHPs_TEST[tit], float(BONDS_TEST[it].split()[5])))
        ftw.write("%d %d %10f %10f\n" % (tit, it, ICOHPs_TEST[tit], float(BONDS_TEST[it].split()[5])))
        dtpx.append(ICOHPs_TEST[tit])
        if higher < ICOHPs_TEST[tit]: higher = ICOHPs_TEST[tit]
        dtpy.append(float(BONDS_TEST[it].split()[5]))
        if higher < float(BONDS_TEST[it].split()[5]): higher = float(BONDS_TEST[it].split()[5])
    #else:
    #    print tit, it, ICOHPs_TEST[tit], float(BONDS_TEST[it].split()[5])
    tit += 1

if DIR_TO_WRITE != "":
    dattw.close()
    # writting the root-mean square error ... 
    RMSE = np.sqrt(SUMDIFF/len(ICOHPs_TEST))
    # ... and the variances of training and testing sets
    ftw.write("TRAINING_VAR = %10f\nTRAINING_STD = %10f\nTESTING_VAR = %10f\nTESTING_STD = %10f\nRMSE = %10f\nMAX_ICOHP = %10f\n" 
              % (VAR, STD, np.var(dtpy), np.std(dtpy), RMSE, higher))
    ftw.close()
else:
    # writting the root-mean square error ... 
    RMSE = np.sqrt(SUMDIFF/len(ICOHPs_TEST))
    # ... and the variances of training and testing sets
    print "RMSE = %10f\nTRAINING_VAR = %10f\nTRAINING_STD = %10f\n" % (RMSE, VAR, STD)

# showing and saving the results
if DIR_TO_WRITE != "" and (SHOW_PLOT or SAVE_PLOT):
    fig = plt.figure(num=FILES_PREFIX)
    plt.subplot(1, 1, 1)
    plt.xlim(0,math.ceil(higher*1.2))
    plt.ylim(0,math.ceil(higher*1.2))
    plt.scatter(dtpx, dtpy, color = "k")
    plt.xlabel("predicted -ICOHP (eV)", {"color": "k", "fontsize": 15, "weight" : "bold"})
    plt.ylabel("target -ICOHP (eV)", {"color": "k", "fontsize": 15, "weight" : "bold"})
    plt.title("%s-%s" % 
              (SPECIE_A, SPECIE_B), 
              {"color": "k", "fontsize": 15, "weight" : "bold"})
    if sys.argv[6+iarg] == "SE" or sys.argv[6+iarg] == "OU":
        plt.text(0.15, 0.61, 
                 "training set size: %5d (var = %6.4f, std = %6.4f)\ntesting set size:  %5d (var = %6.4f, std = %6.4f)\n$\zeta$ = %s\n$\\theta$ = %4.2f\n$\gamma$ = %5.3f\nRMSE = %6.4f\nmax(-ICOHP) = %4.2f" % 
                 (KERNEL_DIM, VAR, STD, TEST_SIZE, np.var(dtpy), np.std(dtpy), ZETA, SIGMA, REGP, RMSE, higher), 
                 {"color": "k", "fontsize": 10, "weight" : "bold"}, transform = fig.transFigure)
    else:
        plt.text(0.15, 0.61, 
                 "training set size: %5d (var = %6.4f, std = %6.4f)\ntesting set size:  %5d (var = %6.4f, std = %6.4f)\n$\zeta$ = %4.2f\n$\\theta$ = %4.2f\n$\gamma$ = %5.3f\nRMSE = %6.4f\nmax(-ICOHP) = %4.2f" % 
                 (KERNEL_DIM, VAR, STD, TEST_SIZE, np.var(dtpy), np.std(dtpy), ZETA, SIGMA, REGP, RMSE, higher), 
                 {"color": "k", "fontsize": 10, "weight" : "bold"}, transform = fig.transFigure)
    if SAVE_PLOT: plt.savefig(DIR_TO_WRITE+"/"+FILES_PREFIX+".png")
    if SHOW_PLOT: plt.show()

