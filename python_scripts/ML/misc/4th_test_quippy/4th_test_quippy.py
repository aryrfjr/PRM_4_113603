#!/usr/bin/python

# http://libatoms.github.io/QUIP/Tutorials/Introduction.html
# http://libatoms.github.io/QUIP/Tutorials/quippy-descriptor-tutorial.html#A-many-body-descriptor:-SOAP

import numpy as np
import quippy
from quippy import Atoms
from quippy import descriptors

#
def create_desc(at, argsd):
    # only Al local environments
    desc = descriptors.Descriptor(argsd)
    at.set_cutoff(desc.cutoff())
    at.calc_connect()
    soap = desc.calc(at, grad=False)["descriptor"]
    centers = desc.calc(at, grad=False)["descriptor_index_0based"]
    return soap, centers

#
def write_lobsterin(catis, atoms):
    ftw = open('lobsterin','w')
    ftw.write("COHPstartEnergy -80.0\n")
    ftw.write("COHPendEnergy 40.0\n")
    ftw.write("basisSet Bunge\n")
    ftw.write("basisfunctions Al 3s 3p\n")
    ftw.write("basisfunctions Cu 4s 4p 3d\n")
    ftw.write("basisfunctions Zr 4s 5s 4p 5p 4d\n")
    csymbs = atoms.get_chemical_symbols()
    for k in range(len(catis)):
        indices, offsets = atoms.connect.get_neighbours(catis[k] + 1)
        for i, offset in zip(indices, offsets):
            ftw.write("cohpbetween atom %d atom %d\n" % (catis[k] + 1, i))
    ftw.close()

########################################################################
# the script
########################################################################
# python ~/dev-python/utils/dump2/dump2xyz.py /home/aryjr/UFSCar/MG-NMR/BMG/Zr49Cu49Al2/c/md/lammps/2/2 151000000 cell
# python ~/dev-python/utils/dump2/dump2xyz.py /home/aryjr/UFSCar/MG-NMR/BMG/Zr49Cu49Al2/c/md/lammps/2/2 151000000 coord

lmax = 4 # spherical harmonics basis band limit
nmax = 5 # number of radial basis functions
atl = quippy.AtomsList("Zr49Cu49Al2.xyz")
# cutoff from RDFs reported in [PRL 102, 245501 (2009)]
args = "soap cutoff=3.7 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={40} n_species=3 species_Z={13 29 40} average=True"

# Computing the SOAP vectors
soap, cats = create_desc(atl[len(atl)-1], args)
atoms = atl[len(atl)-1]
write_lobsterin(cats[0], atoms)

