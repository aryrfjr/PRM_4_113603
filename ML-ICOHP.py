#!/usr/bin/python
#

#
# This is a script that reads a set of steps from a 
# LAMMPS .dump file to compute the ML predicted ICOHP 
# values for all detected bonds.
#

# Libraries
import numpy as np
import pickle as pkl
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_lammps_dump
from theo4m.atom.desc.soap import soaplist
from theo4m.atom.interaction import Interaction

##########################################################################
#
# Subroutines
#
##########################################################################

# Loads the database
def load_database():
    DB_DIR = BIGDATA_DIR+DB_NAME
    for iildb in range(len(INTERACTIONS)):
        print ("Loading database for interaction %s:\n" % INTERACTIONS[iildb])
        # Loading the bonds of the training set and their respective SOAPs
        fileobj = open(DB_DIR+"/"+DB_NAME+"_"+INTERACTIONS[iildb]+".bnd")
        BONDS_TRAIN.append(fileobj.readlines())
        fileobj.close()
        ###########################################################################
        # previous version using Python 2 (see on 27/06/2020-(2))
        #ftl = open(DB_DIR+"/"+DB_NAME+"_"+INTERACTIONS[iildb]+"-SOAPS.vec","rb")
        #SOAPS_TRAIN.append(pkl.load(ftl))
        ###########################################################################
        # https://stackoverflow.com/questions/50283123/python-3-pickle-load-from-python-2
        with open(DB_DIR+"/"+DB_NAME+"_"+INTERACTIONS[iildb]+"-SOAPS.vec", "rb") as ftl:
            SOAPS_TRAIN.append(pkl.load(ftl, fix_imports=True, encoding="latin1"))
        ftl.close()

# Computes per-interaction alpha_i (eq. 3 in Supp. Info. of [Nature Comm. 9, 4501 (2018)])
def compute_alpha(id_int):
    print ("Computing alpha_i for interaction %s:" % INTERACTIONS[id_int])
    # Building the kernel matrix
    btint = BONDS_TRAIN[id_int]
    stint = SOAPS_TRAIN[id_int]
    kd = len(btint) # the kernel dimension
    KM = np.zeros((kd,kd)) # the kernel matrix
    for ik in range(kd):
        bi_dist = float(btint[ik].split()[5])
        soapAi = stint[btint[ik].split()[0]+"-"+btint[ik].split()[1]+"-"+btint[ik].split()[2]+"-"+btint[ik].split()[3]]
        soapBi = stint[btint[ik].split()[0]+"-"+btint[ik].split()[1]+"-"+btint[ik].split()[2]+"-"+btint[ik].split()[4]]
        for jk in range(kd):
            bj_dist = float(btint[jk].split()[5])
            soapAj = stint[btint[jk].split()[0]+"-"+btint[jk].split()[1]+"-"+btint[jk].split()[2]+"-"+btint[jk].split()[3]]
            soapBj = stint[btint[jk].split()[0]+"-"+btint[jk].split()[1]+"-"+btint[jk].split()[2]+"-"+btint[jk].split()[4]]
            dpAAij = np.dot(soapAi,soapAj)
            dpBBij = np.dot(soapBi,soapBj)
            dpABij = np.dot(soapAi,soapBj)
            dpBAij = np.dot(soapBi,soapAj)
            KM[ik][jk] = np.exp((-(bi_dist-bj_dist)**2)/(2*SIGMA[id_int]**2))*(0.25*(dpAAij+dpBBij+dpABij+dpBAij))**ZETA[id_int]
    # Now the variance and standard deviation of the -ICOHP values in the training set
    ICOHPs_TRAIN.append([])
    for ik in range(kd):
        ICOHPs_TRAIN[id_int].append(float(btint[ik].split()[6]))
    STD_TRAIN.append(np.std(ICOHPs_TRAIN[id_int]))
    VAR_TRAIN.append(STD_TRAIN[id_int]**2)
    print ("    training set size: %5d (var = %6.4f, std = %6.4f)\n" % (kd, VAR_TRAIN[id_int], STD_TRAIN[id_int]))
    # Next the inverse kernel matrix ... 
    ireg = STD_TRAIN[id_int] * REG_PAR[id_int] * np.identity(kd)
    ikm = np.linalg.inv(KM + ireg)
    # ... and finally the weights (alpha_i)
    ALPHA.append([0.0] * kd)
    for ia in range(kd):
        for ja in range(kd):
            ALPHA[id_int][ia] += ikm.item((ia,ja)) * ICOHPs_TRAIN[id_int][ja]

# Computes the ML -ICOHP for a bond
def get_ml_icohp(bondp, ii):
    # Bond information
    iA = int(bondp.get_symbA()[2:])
    iB = int(bondp.get_symbB()[2:])
    # Selecting the database
    btint = BONDS_TRAIN[ii]
    stint = SOAPS_TRAIN[ii]
    kd = len(btint) # the kernel dimension
    # The bond feature vector components
    bit_dist = bondp.get_distance()
    soapAit = SOAPS_TEST[iA-1]
    soapBit = SOAPS_TEST[iB-1]
    # Predicting the -ICOHP value
    icohpp = 0.0
    for ik in range(kd):
        bi_dist = float(btint[ik].split()[5])
        soapAi = stint[btint[ik].split()[0]+"-"+btint[ik].split()[1]+"-"+btint[ik].split()[2]+"-"+btint[ik].split()[3]]
        soapBi = stint[btint[ik].split()[0]+"-"+btint[ik].split()[1]+"-"+btint[ik].split()[2]+"-"+btint[ik].split()[4]]
        dpAAiti = np.dot(soapAit,soapAi)
        dpBBiti = np.dot(soapBit,soapBi)
        dpABiti = np.dot(soapAit,soapBi)
        dpBAiti = np.dot(soapBit,soapAi)
        # Here I"m using ZETA as suggested by Prof. Gabor on 05/07/2019
        KERNEL = np.exp((-(bit_dist-bi_dist)**2)/(2*SIGMA[ii]**2))*(0.25*(dpAAiti+dpBBiti+dpABiti+dpBAiti))**ZETA[ii]
        # Here I already used ZETA before, as suggested by Prof. Gabor on 05/07/2019
        icohpp += ALPHA[ii][ik]*KERNEL
    return icohpp

##########################################################################
#
# The script
#
##########################################################################

# Constants
BIGDATA_DIR = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/"
REF_SYMBS = ["Zr","Cu","Al"] # the output must be in the order of LAMMPS
# Per-interaction quantities
INTERACTIONS = ["Al-Al", "Cu-Al", "Cu-Cu", "Zr-Al", "Zr-Cu", "Zr-Zr"]
DMIN = [2.1, 2.0, 2.0, 2.5, 2.3, 2.5] # minimal bond distance
DMAX = [3.5, 3.5, 3.5, 4.0, 3.9, 4.2] # maximum bond distance
ALPHA = [] # alpha_i (eq. 3 in Supp. Info. of [Nature Comm. 9, 4501 (2018)])
BONDS_TRAIN = []
SOAPS_TRAIN = []
ICOHPs_TRAIN = []
STD_TRAIN = []
VAR_TRAIN = []
ZETA = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
SIGMA = [0.5, 0.5, 0.5, 0.5, 1.0, 0.5]
REG_PAR = [0.005, 0.005, 0.010, 0.010, 0.040, 0.010]
# Input arguments
DUMP_FILE = sys.argv[1]
DB_NAME = sys.argv[2]
WORK_DIR = sys.argv[3]
ITEMS = np.array(sys.argv[4:]).astype(int).tolist()

# Loading the databases
load_database()

# Computing per-interaction alpha_i (eq. 3 in Supp. Info. of [Nature Comm. 9, 4501 (2018)])
for ii in range(len(INTERACTIONS)):
    compute_alpha(ii)

# Reading the steps
isteps, iatoms = read_lammps_dump(lmpdump_file = DUMP_FILE, 
                                     spc_symbs = REF_SYMBS, 
                                     frac = True, items = ITEMS) # items = [-1] for all items

# And looping over all of them
for iit in range(len(ITEMS)):
    atoms = iatoms[iit]
    ccell = atoms.get_cell()
    cv = [ccell[0][0], ccell[1][1], ccell[2][2]]
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()

    # And now the SOAP vectors
    soapsl = soaplist.SOAPList(iatoms, verb = True)
    soapsl.compute_per_atom(cutoff = "3.75", 
                  l_max = "6", 
                  n_max = "8", 
                  n_Z = "3", 
                  Z = "{13 29 40}", 
                  n_species = "3", 
                  species_Z = "{13 29 40}")
    # Taking only the per-atom soaps and the 
    # respective central atoms indexes (cats) 
    # and AtomsList object (frame) created 
    # with the only set of arguments above. 
    # Only the last step will be used
    SOAPS_TEST, cats = soapsl.get_per_atom_soaps(0,0)

    # Now looping over the atoms to detect the bonds
    idb = 1
    ftw = open(WORK_DIR+"/"+str(ITEMS[iit])+"_"+isteps[iit]+"_ICOHPLIST.lobster","w")
    ftw.write("  COHP#  atomMU  atomNU     distance     translation      ICOHP (at) eF   for spin  1\n")
    COUNT_INT = [0] * len(INTERACTIONS)
    for i in range(len(atoms)):
        for j in range(i, len(atoms)):
            dv = geom[i] - geom[j]
            trans = [0,0,0]
            for ic in range(3): # PBC
                if abs(dv[ic]) > cv[ic]/2:
                    if dv[ic] > 0:
                        dv[ic] -= cv[ic]
                        trans[ic] = -1
                    elif dv[ic] < 0:
                        dv[ic] += cv[ic]
                        trans[ic] = 1
            dist = np.sqrt(dv[0]**2 + dv[1]**2 + dv[2]**2)
            try:
                iint = INTERACTIONS.index(symbs[i]+"-"+symbs[j])
            except ValueError:
                iint = INTERACTIONS.index(symbs[j]+"-"+symbs[i])
            if dist >= DMIN[iint] and dist <= DMAX[iint]:
                COUNT_INT[iint] += 1
                bond = Interaction()
                bond.set_id(idb)
                # see commentaries on 29/06/2020-(4)
                bond.set_symbA(symbs[i]+str(i+1))
                bond.set_symbB(symbs[j]+str(j+1))
                bond.set_distance(dist)
                bond.set_transx(trans[0])
                bond.set_transy(trans[1])
                bond.set_transz(trans[2])
                # Predicting its -ICOHP value
                bond.set_info(get_ml_icohp(bond, iint))
                #bond.set_info(0.0)
                idb += 1
                # And writing its +ICOHP value
                ftw.write("%7d %7s %7s %12.5f %7d %3d %3d %18.5f\n" % (bond.get_id(), 
                      bond.get_symbA(), bond.get_symbB(), bond.get_distance(), 
                      bond.get_transx(), bond.get_transy(), 
                      bond.get_transz(), -bond.get_info()))
    ftw.close()
    # writing some additional information
    ftw = open(WORK_DIR+"/"+str(ITEMS[iit])+"_"+isteps[iit]+"_ICOHPLIST.info","w")
    ftw.write("%d interactions\n\n" % (idb - 1))
    ftw.write("[interaction] [count] [min. dist.] [max. dist.] [zeta] [sigma] [reg. par.]\n")
    for iint in range(len(INTERACTIONS)):
        ftw.write("%10s %10d  %15.2f %15.2f %15.2f %15.2f %18.5f\n" % (INTERACTIONS[iint], 
              COUNT_INT[iint], DMIN[iint], DMAX[iint], 
              ZETA[iint], SIGMA[iint], REG_PAR[iint]))
    ftw.close()

