#!/usr/bin/python

# References:
# http://libatoms.github.io/QUIP/Tutorials/Introduction.html
# http://libatoms.github.io/QUIP/Tutorials/quippy-descriptor-tutorial.html#A-many-body-descriptor:-SOAP

import sys
import numpy as np
import quippy
#from quippy import Atoms
from quippy import descriptors

# create a global SOAP descriptor for a given single-frame .xyz file
def create_desc(at, argsd, debug=False):
    desc = descriptors.Descriptor(argsd)
    at.set_cutoff(desc.cutoff())
    at.calc_connect()
    soap = desc.calc(at, grad=False)["descriptor"]
    return soap

########################################################################
# the script
########################################################################

lmax = 6 # spherical harmonics basis band limit
nmax = 8 # number of radial basis functions

CHEM_COMPOSITION_A = sys.argv[1]
NAT_A = sys.argv[2]
ID_RUN_A = sys.argv[3]
path_A = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/"+CHEM_COMPOSITION_A+"/c/md/lammps/"+NAT_A+"/"+ID_RUN_A+"/"
STEP_A = sys.argv[4]

CHEM_COMPOSITION_B = sys.argv[5]
NAT_B = sys.argv[6]
ID_RUN_B = sys.argv[7]
path_B = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/"+CHEM_COMPOSITION_B+"/c/md/lammps/"+NAT_B+"/"+ID_RUN_B+"/"
STEP_B = sys.argv[8]

# reading the single-frame ([0] index) .xyz files
ATOMS_A = quippy.AtomsList(path_A+"extxyz/"+STEP_A+".xyz")[0]
ATOMS_B = quippy.AtomsList(path_B+"extxyz/"+STEP_B+".xyz")[0]

SOAPS_A = []
SOAPS_B = []
# Zr central atoms
args = "soap cutoff=3.75 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={40} n_species=3 species_Z={13 29 40} average=True"
SOAPS_A.append(create_desc(ATOMS_A, args))
SOAPS_B.append(create_desc(ATOMS_B, args))
# Cu central atoms
args = "soap cutoff=3.0 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={29} n_species=3 species_Z={13 29 40} average=True"
SOAPS_A.append(create_desc(ATOMS_A, args))
SOAPS_B.append(create_desc(ATOMS_B, args))
# Al central atoms
args = "soap cutoff=3.0 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={13} n_species=3 species_Z={13 29 40} average=True"
SOAPS_A.append(create_desc(ATOMS_A, args))
SOAPS_B.append(create_desc(ATOMS_B, args))
# Zr, Cu, and Al central atoms
args = "soap cutoff=3.75 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=3 Z={13 29 40} n_species=3 species_Z={13 29 40} average=True"
SOAPS_A.append(create_desc(ATOMS_A, args))
SOAPS_B.append(create_desc(ATOMS_B, args))

print "%s %10f %10f %10f %10f" % (STEP_B, np.dot(SOAPS_A[0][0],SOAPS_B[0][0]), np.dot(SOAPS_A[1][0],SOAPS_B[1][0]), np.dot(SOAPS_A[2][0],SOAPS_B[2][0]), np.dot(SOAPS_A[3][0],SOAPS_B[3][0]))

