#!/usr/bin/python

#
# Compares clusters in an extended .xyz file using their SOAPS.
# 
# usage:
# 
# python compare_clusters_soap.py <XYZ_FILE> <ID_CLUSTER_A> <ID_CLUSTER_B>
#

# Libraries
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.atom.desc.soap import soaplist
import numpy as np

xyz_file = sys.argv[1]
ida = int(sys.argv[2])
idb = int(sys.argv[3])

soapsl = soaplist.SOAPList(xyz_file, verb = True)
soapsl.compute_per_atom(cutoff = "3.75", 
                  l_max = "6", 
                  n_max = "8", 
                  n_Z = "3", 
                  Z = "{13 29 40}", 
                  n_species = "3", 
                  species_Z = "{13 29 40}")

soapsa, catsa, dummy = soapsl.get_pasoaps(0,ida)
soapsb, catsb, dummy = soapsl.get_pasoaps(0,idb)
for i in range(len(soapsa)):
    print catsa[i][0], catsb[i][0], np.dot(soapsa[i],soapsb[i])

