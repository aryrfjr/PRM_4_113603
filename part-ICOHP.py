#!/usr/bin/python
#

#
# See on 12/08/2019-(9).
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
OUT_DIR = sys.argv[4]
ICOHP_SPLIT = float(sys.argv[5])

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
icohp_atoms = defaultdict(list)
fcs_atoms = defaultdict(list)
for bk, b in bonds.iteritems():
    # filling up a per-atom dictionary with the accumulated ICOHP values for each atom
    icohp_atoms[b.get_symbA()].append(-b.get_info())
    icohp_atoms[b.get_symbB()].append(-b.get_info())
    # filling up a per-atom dictionary with the symbols of the FCS atoms
    fcs_atoms[b.get_symbA()].append(b.get_symbB()[:2])
    fcs_atoms[b.get_symbB()].append(b.get_symbA()[:2])

# Now dividing the groups
groups = {}
groups["G1Zr"] = []
groups["G2Zr"] = []
groups["G1Cu"] = []
groups["G2Cu"] = []
groups["G1Al"] = []
groups["G2Al"] = []
for ka, sicohp in icohp_atoms.iteritems():
    if sum(sicohp) / len(sicohp) <= ICOHP_SPLIT:
        if ka[:2] == "Zr": 
            groups["G1Zr"].append(int(ka[2:])-1)
        elif ka[:2] == "Cu":
            groups["G1Cu"].append(int(ka[2:])-1)
        elif ka[:2] == "Al":
            groups["G1Al"].append(int(ka[2:])-1)
    else:
        if ka[:2] == "Zr": 
            groups["G2Zr"].append(int(ka[2:])-1)
        elif ka[:2] == "Cu":
            groups["G2Cu"].append(int(ka[2:])-1)
        elif ka[:2] == "Al":
            groups["G2Al"].append(int(ka[2:])-1)

# Writing .xyz and .grp files for each group
for kg, grp in groups.iteritems():
    xyz_file = OUT_DIR+"/"+str(ITEM)+"_"+isteps[0]+"_"+kg+".xyz"
    print xyz_file
    ftw = open(xyz_file,"w")
    ftw.write("%d\n" % len(grp))
    # The order of the lattice vector was wrong in a previous version 
    # but not now it is correct. See more at:
    # https://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz
    ftw.write("Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3:pZr:R:1:pCu:R:1:pAl:R:1:aicohp:R:1 pbc=\"T T T\"\n" 
              % (ccell[0,0], ccell[0,1], ccell[0,2], ccell[1,0], ccell[1,1], ccell[1,2], ccell[2,0], ccell[2,1], ccell[2,2]))
    for i in grp:
        fcs_ca = fcs_atoms[symbs[i]+str(i+1)] # the first coordination shell
        # percentual composition of the FCS
        pZr = float(fcs_ca.count("Zr")) / float(len(fcs_ca)) * 100.0
        pCu = float(fcs_ca.count("Cu")) / float(len(fcs_ca)) * 100.0
        pAl = float(fcs_ca.count("Al")) / float(len(fcs_ca)) * 100.0
        sicohp = icohp_atoms[symbs[i]+str(i+1)]
        aicohp = sum(sicohp) / len(sicohp) # the average -ICOHP values
        ftw.write("%s %10f %10f %10f %10f %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2], pZr, pCu, pAl, aicohp))
    ftw.close()

    grp_file = OUT_DIR+"/"+str(ITEM)+"_"+isteps[0]+"_"+kg+".grp"
    print grp_file
    ftw = open(grp_file,"w")
    for i in grp: ftw.write(str(ids[i])+" ")
    ftw.close()

