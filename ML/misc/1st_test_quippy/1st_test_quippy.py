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

'''
# From the quippy source code file: libmatch/soap.py, subroutine Soap2AlchemySoap
# 
def quippy_2_alchemSOAP(rawsoap, lmax, nmax, n_species, species_Z):
    dimas = nmax ** 2 * (lmax + 1) # dimension of the alchemical SOAP vector
    asdict = {} # dictionary with tuples (Z1,Z2)
    ipair = {}
    # initialize the alchemical soap
    for s1 in xrange(n_species):
        for s2 in xrange(n_species):
            asdict[(species_Z[s1], species_Z[s2])] = np.zeros(dimas, float)
            ipair[(species_Z[s1], species_Z[s2])] = 0
    # 
    isoap = 0
    isqrttwo = 1.0 / np.sqrt(2.0) 
    # selpair and revpair are modified and in turn modify soaps because they are all pointing at the same memory block
    for s1 in xrange(n_species):
        for n1 in xrange(nmax): # loop over the number of radial basis functions
            for s2 in xrange(s1 + 1):
                selpair = asdict[(species_Z[s2], species_Z[s1])]
                # we need to reconstruct the spectrum for the inverse species order, that also swaps n1 and n2.
                # This is again only needed to enable alchemical combination of e.g. alpha-beta and beta-alpha. Shit happens
                revpair = asdict[(species_Z[s1], species_Z[s2])]
                isel = ipair[(species_Z[s2], species_Z[s1])]
                for n2 in xrange(nmax if s2 < s1 else n1 + 1):
                    for l in xrange(lmax + 1):
                        print s1, s2, n1, n2, isel, l+(lmax+1)*(n2+nmax*n1), l+(lmax+1)*(n1+nmax*n2)
                        if (s1 != s2):
                            # undo the normalization since we will actually sum over all pairs in all directions!
                            selpair[isel] = rawsoap[isoap] * isqrttwo
                            revpair[l + (lmax + 1) * (n1 + nmax * n2)] = selpair[isel]
                        else:
                            # diagonal species (s1=s2) have only half of the elements.
                            # this is tricky. we need to duplicate diagonal blocks "repairing" these to be full.
                            # this is necessary to enable alchemical similarity matching, where we need to combine
                            # alpha-alpha and alpha-beta environment fingerprints
                            selpair[l + (lmax + 1) * (n2 + nmax * n1)] = rawsoap[isoap]
                            selpair[l + (lmax + 1) * (n1 + nmax * n2)] = rawsoap[isoap]
                        isoap += 1
                        isel += 1
                ipair[(species_Z[s2], species_Z[s1])] = isel
    return asdict
'''

########################################################################
# the script
########################################################################

# a test set of frames starting from the perfect conventional cell 
# of ZrCu2Al followed by exchanges and removal of Zr and Cu atoms
lmax = 4 # spherical harmonics basis band limit
nmax = 4 # number of radial basis functions
at_z2ca = quippy.AtomsList("Zr2CuAl-conv.xyz")
at_zc2a = quippy.AtomsList("ZrCu2Al-conv.xyz")
at_both = quippy.AtomsList("Zr2CuAl_ZrCu2Al-conv.xyz")
at_z2cad = quippy.AtomsList("Zr2CuAl-conv_dist.xyz")
at_zc2ad = quippy.AtomsList("ZrCu2Al-conv_dist.xyz")
at_zc2amd = quippy.AtomsList("ZrCu2Al-conv_mdist.xyz")
argsd = "soap cutoff=3.1 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=1 Z={13} n_species=3 species_Z={13 29 40} average=False"

# One QUIP SOAP spectrum and the respective atoms for each frame
qsoaps = []
cats = []
for iframe in range(len(at_z2ca)):
    dsoap, dcats = create_desc(at_z2ca[iframe], argsd)
    qsoaps.append(dsoap)
    cats.append(dcats)

for iframe in range(len(at_zc2a)):
    dsoap, dcats = create_desc(at_zc2a[iframe], argsd)
    qsoaps.append(dsoap)
    cats.append(dcats)

for iframe in range(len(at_both)):
    dsoap, dcats = create_desc(at_both[iframe], argsd)
    qsoaps.append(dsoap)
    cats.append(dcats)

for iframe in range(len(at_z2cad)):
    dsoap, dcats = create_desc(at_z2cad[iframe], argsd)
    qsoaps.append(dsoap)
    cats.append(dcats)

for iframe in range(len(at_zc2ad)):
    dsoap, dcats = create_desc(at_zc2ad[iframe], argsd)
    qsoaps.append(dsoap)
    cats.append(dcats)

for iframe in range(len(at_zc2amd)):
    dsoap, dcats = create_desc(at_zc2amd[iframe], argsd)
    qsoaps.append(dsoap)
    cats.append(dcats)

'''
print 'Global comparison between frames'
for i in range(6):
    for j in range(i + 1, 6):
        print i, j, np.dot(qsoaps[i][0], qsoaps[j][0])
print ''
'''

print 'Environments in Zr2CuAl-conv.xyz'
for i in range(len(qsoaps[0])):
    for j in range(i + 1, len(qsoaps[0])):
        print i, j, np.dot(qsoaps[0][i], qsoaps[0][j])
print ''

print 'Environments in ZrCu2Al-conv.xyz'
for i in range(len(qsoaps[1])):
    for j in range(i + 1, len(qsoaps[1])):
        print i, j, np.dot(qsoaps[1][i], qsoaps[1][j])
print ''

print 'Environments in the 1st frame of Zr2CuAl_ZrCu2Al-conv.xyz'
for i in range(len(qsoaps[2])):
    for j in range(i + 1, len(qsoaps[2])):
        print i, j, np.dot(qsoaps[2][i], qsoaps[2][j])
print ''

print 'Environments in the 2nd frame of Zr2CuAl_ZrCu2Al-conv.xyz'
for i in range(len(qsoaps[3])):
    for j in range(i + 1, len(qsoaps[3])):
        print i, j, np.dot(qsoaps[3][i], qsoaps[3][j])
print ''

print 'Environments in Zr2CuAl-conv_dist.xyz'
for i in range(len(qsoaps[4])):
    for j in range(i + 1, len(qsoaps[4])):
        print i, j, np.dot(qsoaps[4][i], qsoaps[4][j])
print ''

print 'Environments in ZrCu2Al-conv_dist.xyz'
for i in range(len(qsoaps[5])):
    for j in range(i + 1, len(qsoaps[5])):
        print i, j, np.dot(qsoaps[5][i], qsoaps[5][j])
print ''

print 'The equivalent environments in Zr2CuAl-conv.xyz and in ZrCu2Al-conv.xyz'
for i in range(len(qsoaps[0])):
    print i, np.dot(qsoaps[0][i], qsoaps[1][i])
print ''

print 'The equivalent environments in Zr2CuAl-conv.xyz and in the 1st frame of Zr2CuAl_ZrCu2Al-conv.xyz'
for i in range(len(qsoaps[0])):
    print i, np.dot(qsoaps[0][i], qsoaps[2][i])
print ''

print 'The equivalent environments in Zr2CuAl-conv.xyz and in the 2nd frame of Zr2CuAl_ZrCu2Al-conv.xyz'
for i in range(len(qsoaps[0])):
    print i, np.dot(qsoaps[0][i], qsoaps[3][i])
print ''

print 'The equivalent environments in Zr2CuAl-conv.xyz and in Zr2CuAl-conv_dist.xyz'
for i in range(len(qsoaps[0])):
    print i, np.dot(qsoaps[0][i], qsoaps[4][i])
print ''

print 'The equivalent environments in Zr2CuAl-conv.xyz and in ZrCu2Al-conv_dist.xyz'
for i in range(len(qsoaps[0])):
    print i, np.dot(qsoaps[0][i], qsoaps[5][i])
print ''

print 'The equivalent environments in ZrCu2Al-conv_dist.xyz and in ZrCu2Al-conv_mdist.xyz'
for i in range(len(qsoaps[0])):
    print i, np.dot(qsoaps[5][i], qsoaps[6][i])
print ''

'''
# Now computing the alchemical soaps for each center in each frame
asoaps = []
for iframe in range(len(at_cswp)):
    asoaps.append([])
    tot = 0
    for icat in range(len(qsoaps[iframe])):
        print len(qsoaps[iframe][icat])
        asoaps[iframe].append(quippy_2_alchemSOAP(qsoaps[iframe][icat], lmax, nmax, 3, [13, 29, 40]))
        print asoaps[iframe][icat].keys()
        for key in asoaps[iframe][icat].keys():
            print key
            print asoaps[iframe][icat][key]
'''

