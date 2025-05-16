#!/usr/bin/python

#
# See more about this script on 17/06/2019-(4).
# 
# python check_ICOHPLIST.py <CHEM_COMPOSITION> <ID_RUN> <STEP> <SUB_STEP> <WRITE>
#

# Libraries
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_cohp
from ase.io import read
import numpy as np

#
# Subroutines
#

# writing an individual extended .xyz file
def eval_cluster(acs, offset):
    symbA = symbs[acs[0]-1]
    knf = "KNF: " # keys not found
    for ib in range(1, len(acs)):
        symbB = symbs[acs[ib]-1]
        key = symbA+str(acs[0])+"-"+symbB+str(acs[ib])
        nknf = ""
        try:
            bond = bonds[key]
        except KeyError:
            nknf = key
            try:
                key = symbB+str(acs[ib])+"-"+symbA+str(acs[0])
                bond = bonds[key]
                nknf = ""
            except KeyError:
                nknf = key
        # if one of the two keys have not been found
        if nknf != "": knf = knf + nknf + ", "
    if knf != "KNF: ": print "Problem!!!", knf
    if write:
        ftw = open(jobdir+"cluster-"+symbA+str(acs[0])+".xyz","w")
        ftw.write("%d\n" % len(acs))
        ftw.write("%s\n" % knf)
        tpos = np.dot(offset, cell)
        for ia in range(len(acs)):
            if ia == 0:
                x = pos[acs[ia]-1][0]
                y = pos[acs[ia]-1][1]
                z = pos[acs[ia]-1][2]
            else:
                x = pos[acs[ia]-1][0] + tpos[ia-1][0]
                y = pos[acs[ia]-1][1] + tpos[ia-1][1]
                z = pos[acs[ia]-1][2] + tpos[ia-1][2]
            ftw.write("%s %10f %10f %10f\n" % (symbs[acs[ia]-1], x, y, z))
        ftw.close()

#
# The script
#

# Constants
BD_DIR = "MG-NMR/ML/big-data-full/"
LOCAL_DIR = "/home/aryjr/UFSCar/"

# Input arguments
nominal_comp = sys.argv[1]
id_run = sys.argv[2]
step = sys.argv[3]
sub_step = sys.argv[4]
if len(sys.argv) > 5:
    write = sys.argv[5] == "write"
else:
    write = False

# Auxiliary constants
jobdir = LOCAL_DIR+BD_DIR+nominal_comp+"/c/md/lammps/100/"+id_run+"/"+step+"/"+sub_step+"/"

# Reading ICOHPLIST.lobster
bonds = read_cohp(cohpfile = jobdir+"ICOHPLIST.lobster", task = "ICOHPLIST")

# Reading the .xyz file
atoms = read(jobdir+nominal_comp+".xyz")
cell = atoms.get_cell()
pos = atoms.get_positions()
symbs = atoms.get_chemical_symbols()

# Reading the reference file with the clusters
fileobj = open(jobdir+"lobsterin-quippy")
lines = fileobj.readlines()
fileobj.close()
lastA = -1
for i in range(9, len(lines)):
    refA = int(lines[i].split()[2])
    refB = int(lines[i].split()[4])
    if refA != lastA:
        if i > 9: eval_cluster(acs, offset)
        acs = []
        acs.append(refA)
        offset = []
    acs.append(refB)
    print lines[i]
    offset.append([int(lines[i].split()[6]), int(lines[i].split()[7]), int(lines[i].split()[8])])
    lastA = refA
# writing the last
eval_cluster(acs, offset)

