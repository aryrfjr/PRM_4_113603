#!/usr/bin/python
#

#
# See on 09/09/2019-(1).
#

# Libraries
import numpy as np
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_lammps
from theo4m.atom.desc.soap import soaplist

# Constants
REF_SYMBS = ["Zr","Cu","Al"] # the output must be in the order of LAMMPS

# Input arguments
DUMP_FILE = sys.argv[1]
REF_FILE = sys.argv[2] # ~/UFSCar/MG-NMR/SS-ML/B2-CuZr/B2-CuZr-333.xyz
WORK_DIR = sys.argv[3] # nucleation
ITEMS = np.array(sys.argv[4:]).astype(int).tolist()

# Loading the SOAP vectors of the reference crystal
# And now the SOAP vectors
soapsl = soaplist.SOAPList(REF_FILE, verb = True)
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
SOAPS_CRYST, cats, frame = soapsl.get_pasoaps(0,0)

# From ~/UFSCar/MG-NMR/SS-ML/B2-CuZr/B2-CuZr-333.xyz
# Taking arbitrarily two atoms at the center of the cell
#                      ID Sym Atm.Num  X/Angstrom     Y/Angstrom     Z/Angstrom 
# ----------------------------------------------------------------------------------
# Selected Atom No.1:  28  Zr  40     +4.893000000   +4.893000000   +4.893000000
# Selected Atom No.2:  27  Cu  29     +3.262000000   +3.262000000   +3.262000000
# ----------------------------------------------------------------------------------
#SOAPS_Zr_CRYST = SOAPS_CRYST[27]
#SOAPS_Cu_CRYST = SOAPS_CRYST[26]

# From ~/UFSCar/MG-NMR/SS-ML/ZrCu2Al/ZrCu2Al-222.xyz
# Taking arbitrarily two atoms at the center of the cell
#                      ID Sym Atm.Num  X/Angstrom     Y/Angstrom     Z/Angstrom 
# ----------------------------------------------------------------------------------
# Selected Atom No.1:  113 Al  13     +6.216300000   +6.216300000   +6.216300000
# Selected Atom No.2:  95  Zr  40     +6.216300000   +3.108150000   +6.216300000
# Selected Atom No.3:  72  Cu  29     +7.770375000   +4.662225000   +4.662225000
# Selected Atom No.4: 
# ----------------------------------------------------------------------------------
SOAPS_Zr_CRYST = SOAPS_CRYST[94]
SOAPS_Cu_CRYST = SOAPS_CRYST[71]
SOAPS_Al_CRYST = SOAPS_CRYST[112]

# Just a test
#isp = 26
#indices, offsets = frame.connect.get_neighbours(cats[isp][0] + 1)
#print isp + 1
#for i, offset in zip(indices, offsets):
#    print i
#sys.exit(0)

# Now the glassy structure
# Reading the steps
isteps, iatoms, iccell = read_lammps(lmpoutput = DUMP_FILE, 
                    spc_symbs = REF_SYMBS, 
                    frac = True, items = ITEMS) # frac = True because the .dump files are not in fractional cooridnates
# And looping over all of them
for iit in range(len(ITEMS)):
    atoms = iatoms[iit]
    ccell = iccell[iit]
    cv = [ccell[0][0], ccell[1][1], ccell[2][2]]
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()

    # Writing a .xyz file to be used by libatoms
    xyz_file = WORK_DIR+"/"+str(ITEMS[iit])+"_"+isteps[iit]+".xyz"
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
    SOAPS_GLASS, cats, frame = soapsl.get_pasoaps(0,0)

    # Looping over the atoms of the glass
    count = 0
    dotref = 0.94
    print "Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"" % (ccell[0,0], ccell[0,1], ccell[0,2], ccell[1,0], ccell[1,1], ccell[1,2], ccell[2,0], ccell[2,1], ccell[2,2])
    for isp in range(len(SOAPS_GLASS)):
        if np.dot(SOAPS_Zr_CRYST,SOAPS_GLASS[isp]) > dotref or np.dot(SOAPS_Cu_CRYST,SOAPS_GLASS[isp]) > dotref or np.dot(SOAPS_Al_CRYST,SOAPS_GLASS[isp]) > dotref:
            count += 1
            print isp, cats[isp], np.dot(SOAPS_Zr_CRYST,SOAPS_GLASS[isp]), np.dot(SOAPS_Cu_CRYST,SOAPS_GLASS[isp]), np.dot(SOAPS_Al_CRYST,SOAPS_GLASS[isp])
            print '%s %16.8f %16.8f %16.8f' % (symbs[isp], geom[isp][0], geom[isp][1], geom[isp][2])
            indices, offsets = frame.connect.get_neighbours(cats[isp][0] + 1)
            for i, offset in zip(indices, offsets):
                print '%s %16.8f %16.8f %16.8f' % (symbs[i-1], geom[i-1][0], geom[i-1][1], geom[i-1][2])
                count += 1
    print count

