#!/usr/bin/python

#
# python mix_SSDBs.py <NDBS> [<NCDBS_1> <NBNDS_1> ...] <MIX_DB_NAME> <SPECIE_A> <SPECIE_B>
#

import os
import sys
import pickle as pkl
from random import shuffle

# Constants
BIGDATA_DIR = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full"

# Reading arguments
NDBS = int(sys.argv[1]) # number of databases
NBNDS = [0] * NDBS # number of bonds from each DB
NCDBS = [""] * NDBS # nominal compositions of each DB
idb = 0
for i in range(2,NDBS*2+1,2):
    NCDBS[idb] = sys.argv[i]
    NBNDS[idb] = int(sys.argv[i+1])
    idb += 1
MIX_DB_NAME = sys.argv[NDBS*2+2]
SPECIE_A = sys.argv[NDBS*2+3]
SPECIE_B = sys.argv[NDBS*2+4]

if not os.path.exists(BIGDATA_DIR+"/"+MIX_DB_NAME):
    os.mkdir(BIGDATA_DIR+"/"+MIX_DB_NAME)

# Now loading the databases
NEW_BNDS = []
NEW_SOAPS = {}
for i in range(NDBS):
    DB_DIR = BIGDATA_DIR+"/"+NCDBS[i]+"-PBSSDB/"
    fileobj = open(DB_DIR+SPECIE_A+"-"+SPECIE_B+".bnd")
    ALL_BNDS = fileobj.readlines()
    fileobj.close()
    ftl = open(DB_DIR+SPECIE_A+"-"+SPECIE_B+"-SOAPS.vec","rb")
    ALL_SOAPS = pkl.load(ftl)
    ftl.close()
    # the original structure of the .bnd files is: <ID_RUN> <SUB_STEP> <ID_ATOM_A> <ID_ATOM_B> <DISTANCE> <-ICOHP>
    # and the keys in ALL_SOAPS are: "<ID_RUN>-<SUB_STEP>-<ID_ATOM>"
    shuffle(ALL_BNDS)
    for j in range(NBNDS[i]):
        bnd = ALL_BNDS[j]
        sb = bnd.split()
        NEW_BNDS.append(str(i)+" "+sb[0]+" "+sb[1]+" "+sb[2]+" "+sb[3]+" "+sb[4]+" "+sb[5])
        soap = ALL_SOAPS[bnd.split()[0]+"-"+bnd.split()[1]+"-"+bnd.split()[2]]
        NEW_SOAPS[str(i)+"-"+bnd.split()[0]+"-"+bnd.split()[1]+"-"+bnd.split()[2]] = soap
        soap = ALL_SOAPS[bnd.split()[0]+"-"+bnd.split()[1]+"-"+bnd.split()[3]]
        NEW_SOAPS[str(i)+"-"+bnd.split()[0]+"-"+bnd.split()[1]+"-"+bnd.split()[3]] = soap

# And finally writing the new database
DB_DIR = BIGDATA_DIR+"/"+MIX_DB_NAME
ftw = open(DB_DIR+"/"+MIX_DB_NAME+"_"+SPECIE_A+"-"+SPECIE_B+".info","w")
ftw.write("%s %s-%s\n" % (MIX_DB_NAME, SPECIE_A, SPECIE_B))
for i in range(NDBS):
    ftw.write("%s %d\n" % (NCDBS[i], NBNDS[i]))
ftw.close()
ftw = open(DB_DIR+"/"+MIX_DB_NAME+"_"+SPECIE_A+"-"+SPECIE_B+".bnd","w")
for i in range(len(NEW_BNDS)):
    ftw.write("%s\n" % NEW_BNDS[i])
ftw.close()
ftd = open(DB_DIR+"/"+MIX_DB_NAME+"_"+SPECIE_A+"-"+SPECIE_B+"-SOAPS.vec","wb")
pkl.dump(NEW_SOAPS, ftd)
ftd.close()

