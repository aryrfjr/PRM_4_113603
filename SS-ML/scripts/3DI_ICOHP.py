#!/usr/bin/python
#

#
# This is a script that simply interpolates the ICOHP values
# localized between two atoms using true values computed with 
# LOBSTER for the 100-atoms cells.
# 

# Libraries
import sys
import numpy as np
from collections import defaultdict
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_cohp
from ase.io import read

#
# The script
#

# Input arguments
FATOMS = "no"
CHEM_COMPOSITION = sys.argv[1]
WTD = sys.argv[2]
ITEM = sys.argv[3]
STEP = sys.argv[4]
if WTD == "atoms": FATOMS = sys.argv[5]
if FATOMS == "yes": GTLT = sys.argv[6]
if FATOMS == "yes": AICOHP_REF = float(sys.argv[7])

# Constants
INP_DIR = "/home/aryjr/UFSCar/MG-NMR/SS-ML/"+CHEM_COMPOSITION+"/c/md/lammps/ultimate/2/icohps_06022020/"
REF_SYMBS = ["Zr","Cu","Al"] # the output must be in the order of LAMMPS

# Reading ICOHPLIST.lobster
bonds = read_cohp(cohpfile = INP_DIR+ITEM+"_"+STEP+"_ICOHPLIST.lobster", task = "ICOHPLIST")
# Reading the .xyz file
atoms = read(INP_DIR+ITEM+"_"+STEP+".xyz")
pos = atoms.get_positions()
symbs = atoms.get_chemical_symbols()
cell = atoms.get_cell()
# Printing the head
print "COUNT_YOURSELF"
print "Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"" % (cell[0,0], cell[0,1], cell[0,2], cell[1,0], cell[1,1], cell[1,2], cell[2,0], cell[2,1], cell[2,2])

if WTD == "bonds": # write bonds with their respective -ICOHP values as H atoms
    bks = bonds.keys()
    for ib in range(len(bks)):
        bond = bonds[bks[ib]]
        pA = [0.0]*3
        # here using the LOBSTER's information about PBC
        pA[0] = pos[int(bks[ib].split("-")[0][2:])-1][0] + cell[0][0] * bond.get_transx()
        pA[1] = pos[int(bks[ib].split("-")[0][2:])-1][1] + cell[1][1] * bond.get_transy()
        pA[2] = pos[int(bks[ib].split("-")[0][2:])-1][2] + cell[2][2] * bond.get_transz()
        pB = pos[int(bks[ib].split("-")[1][2:])-1]
        dvec = pB - pA # the distance vector
        pm = pA + dvec / 2 # half of it is a point between the two bonded atoms
        # now checking whether the point is inside the box
        for ipbc in 0,1,2:
            if pm[ipbc] > cell[ipbc][ipbc]: pm[ipbc] -= cell[ipbc][ipbc]
            if pm[ipbc] < 0.0: pm[ipbc] += cell[ipbc][ipbc]
        if pm[2] > 0.0 and pm[2] <= 2.5:
            print "H %10f %10f %10f %10f" % (pm[0], pm[1], pm[2], -bond.get_info())
elif WTD == "atoms": # write atoms with their respective <-ICOHP> values
    # filling up a per-atom dictionary with the accumulated ICOHP values for each atom
    icohp_atoms = defaultdict(list)
    for bk, b in bonds.iteritems():
        icohp_atoms[b.get_symbA()].append(-b.get_info())
        icohp_atoms[b.get_symbB()].append(-b.get_info())
    # writing the atoms
    if FATOMS == "no":
        for i in range(len(atoms)):
            if pos[i][2] > 0.0 and pos[i][2] <= 2.5:
                icohps = icohp_atoms[symbs[i]+str(i+1)]
                aicohp = sum(icohps) / len(icohps) # the average -ICOHP values
                print "%s %10f %10f %10f %10f" % (symbs[i], pos[i][0], pos[i][1], pos[i][2], aicohp)
    else: # shift atoms with different <-ICOHP>
        for i in range(len(atoms)):
            if pos[i][2] > 0.0 and pos[i][2] <= 2.5:
                icohps = icohp_atoms[symbs[i]+str(i+1)]
                aicohp = sum(icohps) / len(icohps) # the average -ICOHP values
                if GTLT == "lt" and aicohp <= AICOHP_REF:
                    print "%s %10f %10f %10f %10f" % (symbs[i], pos[i][0], pos[i][1], pos[i][2], aicohp)
                elif GTLT == "gt" and aicohp >= AICOHP_REF:
                    print "%s %10f %10f %10f %10f" % (symbs[i], pos[i][0], pos[i][1], pos[i][2], aicohp)

"""
#R = []
#ICOHP = []
...
#R.append(pm)
#ICOHP.append(-bond.get_info())
# writing the atoms
#symbs = atoms.get_chemical_symbols()
#for i in range(len(atoms)):
#    print "%s %10f %10f %10f" % (symbs[i], pos[i][0], pos[i][1], pos[i][2])
"""

"""
# trying https://stackoverflow.com/questions/9419451/3d-contour-plot-from-data-using-mayavi-python
from scipy.interpolate import griddata
import numpy as np
from mayavi.mlab import *
# trying https://stackoverflow.com/questions/9419451/3d-contour-plot-from-data-using-mayavi-python
# Create some test data, 3D gaussian, 200 points
dx, pts = 12, 100j
# Create the grid to interpolate on
X,Y,Z = np.mgrid[0.0:dx:pts, 0.0:dx:pts, 0.0:dx:pts]
# Interpolate the data
F = griddata(R, ICOHP, (X,Y,Z),fill_value=0.0)
# Plot
contour3d(F,contours=1,transparent=True)
"""

