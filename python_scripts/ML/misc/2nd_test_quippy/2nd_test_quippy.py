#!/usr/bin/python

# http://libatoms.github.io/QUIP/Tutorials/Introduction.html
# http://libatoms.github.io/QUIP/Tutorials/quippy-descriptor-tutorial.html#A-many-body-descriptor:-SOAP

import numpy as np
import quippy
from quippy import Atoms
from quippy import descriptors

#
def create_desc(at, argsd, debug=False):
    # only Al local environments
    desc = descriptors.Descriptor(argsd)
    at.set_cutoff(desc.cutoff())
    at.calc_connect()
    soap = desc.calc(at, grad=False)["descriptor"]
    centers = desc.calc(at, grad=False)["descriptor_index_0based"]
    return soap, centers

########################################################################
# the script
########################################################################

# a test set of frames starting from the perfect conventional cell 
# of ZrCu2Al followed by exchanges and removal of Zr and Cu atoms
lmax = 4 # spherical harmonics basis band limit
nmax = 5 # number of radial basis functions
ats = []
ats.append(quippy.AtomsList("Zr49Cu49Al2-big.xyz"))
ats.append(quippy.AtomsList("Zr48Cu48Al4-big.xyz"))
ats.append(quippy.AtomsList("Zr47Cu47Al6-big.xyz"))
ats.append(quippy.AtomsList("Zr46Cu46Al8-big.xyz"))
ats.append(quippy.AtomsList("Zr45Cu45Al10-big.xyz"))
ats.append(quippy.AtomsList("Zr44Cu44Al12-big.xyz"))
ats.append(quippy.AtomsList("Zr43Cu43Al14-big.xyz"))
'''
ats.append(quippy.AtomsList("Zr49Cu49Al2.xyz"))
ats.append(quippy.AtomsList("Zr48Cu48Al4.xyz"))
ats.append(quippy.AtomsList("Zr47Cu47Al6.xyz"))
ats.append(quippy.AtomsList("Zr46Cu46Al8.xyz"))
ats.append(quippy.AtomsList("Zr45Cu45Al10.xyz"))
ats.append(quippy.AtomsList("Zr44Cu44Al12.xyz"))
ats.append(quippy.AtomsList("Zr43Cu43Al14.xyz"))
ats.append(quippy.AtomsList("ZrCu2Al-conv.xyz"))
ats.append(quippy.AtomsList("Zr2CuAl-conv.xyz"))
'''
argsz = "soap cutoff=3.7 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={40} n_species=3 species_Z={13 29 40} average=True"
argsc = "soap cutoff=3.7 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={29} n_species=3 species_Z={13 29 40} average=True"
argsa = "soap cutoff=3.7 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={13} n_species=3 species_Z={13 29 40} average=True"

# Computing the SOAP vectors
fqsoaps = [] # SOAP vectors
fcats = [] # central atoms
for iat in range(len(ats)):
    qsoaps = [] # SOAP vectors
    cats = [] # central atoms
    for iframe in range(len(ats[iat])):
        dsoap, dcats = create_desc(ats[iat][iframe], argsa)
        qsoaps.append(dsoap)
        cats.append(dcats)
    fqsoaps.append(qsoaps)
    fcats.append(cats)

#for iat in range(len(ats)):
#    for iframe in range(len(ats[iat])):
#        for icati in range(len(fcats[iat][iframe])):
#            cati = fcats[iat][iframe][icati][0]
#            atoms = ats[iat][iframe]
#            csymbs = atoms.get_chemical_symbols()
#            nps = atoms.positions[cati]
#            print iat, iframe, cati, csymbs[cati], nps[0], nps[1], nps[2]
#        print ""
#    print "========================================================================"

#cati = fcats[0][0][0][0]
#atoms = ats[0][0]
#csymbs = atoms.get_chemical_symbols()
#nps = atoms.positions[cati]
#print csymbs[cati], nps[0], nps[1], nps[2]
#indices, offsets = atoms.connect.get_neighbours(cati+1)
#for i, offset in zip(indices, offsets):
#   nps = atoms.positions[i-1] + np.dot(offset, atoms.get_cell())
#   print csymbs[i-1], nps[0], nps[1], nps[2]

print 'Global comparison among different compositions'
for iat in range(len(ats)):
    print np.dot(fqsoaps[0][len(fqsoaps[0])-1][0], fqsoaps[iat][len(fqsoaps[iat])-1][0])
print ''

#print 'Single environment in a same composition'
#for iqs in range(1, len(fqsoaps[0])):
#    print np.dot(fqsoaps[0][0][0], fqsoaps[0][iqs][0])
#print ''

