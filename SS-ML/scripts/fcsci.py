#!/usr/bin/python
#

#
# This is a script that reads a step from a LAMMPS 
# .dump file to compute their FCSI indexes (see on 
# 19/09/2019).
#

# Libraries
from collections import defaultdict
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_cohp
from theo4m.io import read_lammps

##########################################################################
#
# Subroutines
#
##########################################################################

##########################################################################
#
# The script
#
##########################################################################

# Constants
REF_SYMBS = ["Zr","Cu","Al"] # the output must be in the order of LAMMPS

# Input arguments
DUMP_FILE = sys.argv[1]
ITEM = int(sys.argv[2])
ICOHP_DIR = sys.argv[3]

# Reading the only step (usually the first)
isteps, iatoms, iccell = read_lammps(lmpoutput = DUMP_FILE, 
            spc_symbs = REF_SYMBS, 
            frac = True, items = [ITEM]) # frac = True because the .dump files are not in fractional cooridnates
atoms = iatoms[0]
ccell = iccell[0]
cv = [ccell[0][0], ccell[1][1], ccell[2][2]]
symbs = atoms.get_chemical_symbols()
geom = atoms.get_positions()
ids = atoms.get_tags()

# Loading the file with the bonds
icohp_ref_file = ICOHP_DIR+"/"+str(ITEM)+"_"+isteps[0]+"_ICOHPLIST.lobster"
print icohp_ref_file
bonds = read_cohp(cohpfile = icohp_ref_file, task = "ICOHPLIST")
fcs_atoms = defaultdict(list)
for bk, b in bonds.iteritems():
    # filling up a per-atom dictionary with the symbols of the FCS atoms
    fcs_atoms[b.get_symbA()].append(b.get_symbB()[:2])
    fcs_atoms[b.get_symbB()].append(b.get_symbA()[:2])

# Now setting the per-atom FCSI indexes
for k_ca, fcs_ca in fcs_atoms.iteritems():
    pZr = float(fcs_ca.count("Zr")) / float(len(fcs_ca)) * 100.0
    pCu = float(fcs_ca.count("Cu")) / float(len(fcs_ca)) * 100.0
    pAl = float(fcs_ca.count("Al")) / float(len(fcs_ca)) * 100.0
    if pZr > 70.0: print "Zr", k_ca, fcs_ca, pZr, pCu, pAl
    if pCu > 70.0: print "Cu", k_ca, fcs_ca, pZr, pCu, pAl
    if pAl > 70.0: print "Al", k_ca, fcs_ca, pZr, pCu, pAl

