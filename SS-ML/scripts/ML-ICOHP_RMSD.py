#!/usr/bin/python
#

#
# This is a script that reads a LAMMPS .dump file 
# and writes a corresponding file with ML predicted 
# ICOHP values for all detected bonds. The user can 
# define which steps will be read as well as the 
# cutoff radius which in principle has to be the same 
# used in the DFT ICOHP simulations. This code also 
# reports the RMSD values for the last step based on 
# its ICOHPLIST.lobster file.
#

# Libraries
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import math
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_cohp
from theo4m.io import read_lammps
from theo4m.atom.desc.soap import soaplist
from theo4m.atom.interaction import Interaction
from ase.io import read

##########################################################################
#
# Subroutines
#
##########################################################################

# Loads the database
def load_database():
    DB_DIR = BIGDATA_DIR+DB_NAME
    for iildb in range(len(INTERACTIONS)):
        print "Loading database for interaction %s:\n" % INTERACTIONS[iildb]
        # Loading the bonds of the training set and their respective SOAPs
        fileobj = open(DB_DIR+"/"+DB_NAME+"_"+INTERACTIONS[iildb]+".bnd")
        BONDS_TRAIN.append(fileobj.readlines())
        fileobj.close()
        ftl = open(DB_DIR+"/"+DB_NAME+"_"+INTERACTIONS[iildb]+"-SOAPS.vec","rb")
        SOAPS_TRAIN.append(pkl.load(ftl))
        ftl.close()

# Computes per-interaction alpha_i (eq. 3 in Supp. Info. of [Nature Comm. 9, 4501 (2018)])
def compute_alpha(id_int):
    print "Computing alpha_i for interaction %s:" % INTERACTIONS[id_int]
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
            if REG_TYPE == "SE":
                KM[ik][jk] = np.exp((-(bi_dist-bj_dist)**2)/(2*SIGMA[id_int]**2)) # Squared Exponential (SE)
            elif REG_TYPE == "OU":
                KM[ik][jk] = np.exp(-abs(bi_dist-bj_dist)/SIGMA[id_int]) # Ornstein-Uhlenbeck (OU)
            elif REG_TYPE == "SOAP":
                KM[ik][jk] = np.exp((-(bi_dist-bj_dist)**2)/(2*SIGMA[id_int]**2))*(0.25*(dpAAij+dpBBij+dpABij+dpBAij))**ZETA[id_int]
    # Now the variance and standard deviation of the -ICOHP values in the training set
    ICOHPs_TRAIN.append([])
    for ik in range(kd):
        ICOHPs_TRAIN[id_int].append(float(btint[ik].split()[6]))
    STD_TRAIN.append(np.std(ICOHPs_TRAIN[id_int]))
    VAR_TRAIN.append(STD_TRAIN[id_int]**2)
    print "    training set size: %5d (var = %6.4f, std = %6.4f)\n" % (kd, VAR_TRAIN[id_int], STD_TRAIN[id_int])
    # Next the inverse kernel matrix ... 
    ireg = STD_TRAIN[id_int] * REG_PAR[id_int] * np.identity(kd)
    ikm = np.linalg.inv(KM + ireg)
    # ... and finally the weights (alpha_i)
    ALPHA.append([0.0] * kd)
    for ia in range(kd):
        for ja in range(kd):
            ALPHA[id_int][ia] += ikm.item((ia,ja)) * ICOHPs_TRAIN[id_int][ja]

# Computes the ML -ICOHP for a bond
def get_ml_icohp(bondp):
    # Bond information
    iA = int(bondp.get_symbA()[2:])
    sA = bondp.get_symbA()[:2]
    iB = int(bondp.get_symbB()[2:])
    sB = bondp.get_symbB()[:2]
    # Selecting the database
    try:
        ii = INTERACTIONS.index(sA+"-"+sB)
    except ValueError:
        ii = INTERACTIONS.index(sB+"-"+sA)
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
        if REG_TYPE == "SE":
            KERNEL = np.exp((-(bit_dist-bi_dist)**2)/(2*SIGMA[ii]**2)) # Squared Exponential (SE)
        elif REG_TYPE == "OU":
            KERNEL = np.exp(-abs(bit_dist-bi_dist)/SIGMA[ii]) # Ornstein-Uhlenbeck (OU)
        elif REG_TYPE == "SOAP":
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
INTERACTIONS = ["Al-Al", "Cu-Al", "Cu-Cu", "Zr-Al", "Zr-Cu", "Zr-Zr"]
TASK_COMPUTE_RMSD = "rmsd" # RMSD values for the last step
TASK_COMPUTE_ICOHP = "icohp" # writes a file with ML predicted ICOHP values (MPI???)
DMIN = 1.80 # minimal bond distance
DMAX = 4.20 # maximum bond distance
# Per-interaction quantities
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
CHEM_COMPOSITION = sys.argv[1]
ID_RUN = sys.argv[2]
TASK = sys.argv[3]
DB_NAME = sys.argv[4]
SUB_STEP = sys.argv[5]
if len(sys.argv) > 6:
    REG_TYPE = sys.argv[6]
else:
    REG_TYPE = "SOAP"

if TASK == "rmsd": # just to perform a test based on 100-atoms CMDD cells
    # Loading the databases
    load_database()

    # Computing per-interaction alpha_i (eq. 3 in Supp. Info. of [Nature Comm. 9, 4501 (2018)])
    for ii in range(len(INTERACTIONS)):
        compute_alpha(ii)

    if SUB_STEP == "0":
        # Reading the last step (items = [-2])
        dfile = BIGDATA_DIR+CHEM_COMPOSITION+"/c/md/lammps/100/"+ID_RUN+"/zca-th300.dump"
        rsteps, ratoms, rccell = read_lammps(lmpoutput = dfile, 
                                     spc_symbs = REF_SYMBS, 
                                     frac = False, items = [-2])
        atoms = ratoms[0]
        ccell = rccell[0]
    else:
        # Reading the first structure in the Extended .xyz file
        xyz_file = BIGDATA_DIR+CHEM_COMPOSITION+"/c/md/lammps/100/"+ID_RUN+"/2000/"+SUB_STEP+"/"+CHEM_COMPOSITION+".xyz"
        atoms = read(xyz_file)
        ccell = atoms.get_cell()

    cv = [ccell[0][0], ccell[1][1], ccell[2][2]]
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()

    # Computing the atomic SOAPs
    if SUB_STEP == "0":
        # Writing a .xyz file to be used by libatoms
        xyz_file = BIGDATA_DIR+CHEM_COMPOSITION+"/c/md/lammps/100/"+ID_RUN+"/last.xyz"
        ftw = open(xyz_file,"w")
        ftw.write("%d\n" % len(atoms))
        # The order of the lattice vector was wrong in a previous version 
        # but not now it is correct. See more at:
        # https://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz
        ftw.write("Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" 
                  % (ccell[0,0], ccell[0,1], ccell[0,2], ccell[1,0], ccell[1,1], ccell[1,2], ccell[2,0], ccell[2,1], ccell[2,2]))
        for i in range(len(atoms)):
            ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
        ftw.close()

    # And now the SOAP vectors
    soapsl = soaplist.SOAPList(xyz_file, verb = True)
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
    SOAPS_TEST, cats, frame = soapsl.get_pasoaps(0,0)

    # Now looping over the atoms to detect the bonds
    idb = 1
    bonds = {}
    print "  COHP#  atomMU  atomNU     distance     translation      ICOHP (at) eF   for spin  1"
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
            if dist >= DMIN and dist <= DMAX:
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
                bond.set_info(get_ml_icohp(bond))
                bonds[bond.get_symbA()+"-"+bond.get_symbB()] = bond
                print "%7d %7s %7s %12.5f %7d %3d %3d %18.5f" % (bond.get_id(), 
                                bond.get_symbA(), bond.get_symbB(), bond.get_distance(), 
                                bond.get_transx(), bond.get_transy(), bond.get_transz(), -bond.get_info())
                idb += 1

    # Loading the bonds of the testing set to compute RMSD
    icohp_ref_file = BIGDATA_DIR+CHEM_COMPOSITION+"/c/md/lammps/100/"+ID_RUN+"/2000/"+SUB_STEP+"/ICOHPLIST.lobster"
    print icohp_ref_file
    bonds_ref = read_cohp(cohpfile = icohp_ref_file, task = "ICOHPLIST")
    if len(bonds) == len(bonds_ref):
        summdiff = [0.0] * len(INTERACTIONS)
        nbnds = [0] * len(INTERACTIONS)
        higher = [0.0] * len(INTERACTIONS)
        dtpx = []#] * len(INTERACTIONS)
        dtpy = []#] * len(INTERACTIONS)
        for ibk in range(len(bonds)):
            bk = bonds.keys()[ibk]
            b = bonds[bk]
            br = bonds_ref[bk]
            # Selecting the database
            sA = b.get_symbA()[:2]
            sB = b.get_symbB()[:2]
            try:
                ii = INTERACTIONS.index(sA+"-"+sB)
            except ValueError:
                ii = INTERACTIONS.index(sB+"-"+sA)
            if higher[ii] < b.get_info(): higher[ii] = b.get_info()
            if ii == 5: dtpx.append(b.get_info())
            if higher[ii] < -br.get_info(): higher[ii] = -br.get_info()
            if ii == 5: dtpy.append(-br.get_info())
            summdiff[ii] += (b.get_info() + br.get_info())**2 # the signal is inverted in the ICOHPLIST.lobster file
            nbnds[ii] += 1
        for ii in range(len(INTERACTIONS)):
            if nbnds[ii] == 0:
                print "\nThere is no %s interaction" % INTERACTIONS[ii]
            elif nbnds[ii] == 1:
                print "\nThere is one single %s interaction and its squared difference is %10f" % (INTERACTIONS[ii], summdiff[ii])
            else:
                rmsd = np.sqrt(summdiff[ii]/nbnds[ii])
                print "\nRMSD for %d %s interactions is %10f" % (nbnds[ii], INTERACTIONS[ii], rmsd)
                if ii == 5:
                    fig = plt.figure(num=INTERACTIONS[ii])
                    plt.subplot(1, 1, 1)
                    plt.xlim(0,4.0)
                    plt.ylim(0,4.0)
                    plt.scatter(dtpx, dtpy, color = "k")
                    plt.xlabel("predicted -ICOHP (eV)", {"color": "k", "fontsize": 15, "weight" : "bold"})
                    plt.ylabel("target -ICOHP (eV)", {"color": "k", "fontsize": 15, "weight" : "bold"})
                    plt.title("%s" % 
                              (INTERACTIONS[ii]), 
                              {"color": "k", "fontsize": 15, "weight" : "bold"})
                    plt.text(0.15, 0.67, 
                             "testing set size:  %5d \n$\zeta$ = %4.2f\n$\sigma$ = %4.2f\nreg. par. = %5.3f\nRMSE = %6.4f" % 
                             (nbnds[ii], ZETA[ii], SIGMA[ii], REG_PAR[ii], rmsd), 
                             {"color": "k", "fontsize": 10, "weight" : "bold"}, transform = fig.transFigure)
                    plt.savefig(INTERACTIONS[ii]+".png")
                    plt.show()
    else:
        "Ops! The number of bonds should be the same!!!"

