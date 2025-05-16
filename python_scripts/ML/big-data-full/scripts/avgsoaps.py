#!/usr/bin/python

import sys
sys.path.append("/home/aryjr/dev-python")
import numpy as np
from theo4m.atom.desc.soap import soaplist

# The script
lmax = 6 # spherical harmonics basis band limit
nmax = 8 # number of radial basis functions

soaps_49492_1 = soaplist.SOAPList("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr49Cu49Al2/c/md/lammps/100/1/extxyzf.xyz", verb = True)
# Zr, Cu, and Al central atoms
soaps_49492_1.compute_average(cutoff = "3.75", 
                      l_max = "6", 
                      n_max = "8", 
                      n_Z = "1", 
                      Z = "{13}", 
                      n_species = "3", 
                      species_Z = "{13 29 40}")

soaps_49492_2 = soaplist.SOAPList("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr49Cu49Al2/c/md/lammps/100/2/extxyzf.xyz", verb = True)
# Zr, Cu, and Al central atoms
soaps_49492_2.compute_average(cutoff = "3.75", 
                      l_max = "6", 
                      n_max = "8", 
                      n_Z = "1", 
                      Z = "{13}", 
                      n_species = "3", 
                      species_Z = "{13 29 40}")

soaps_444412_1 = soaplist.SOAPList("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr44Cu44Al12/c/md/lammps/100/1/extxyzf.xyz", verb = True)
# Zr, Cu, and Al central atoms
soaps_444412_1.compute_average(cutoff = "3.75", 
                      l_max = "6", 
                      n_max = "8", 
                      n_Z = "1", 
                      Z = "{13}", 
                      n_species = "3", 
                      species_Z = "{13 29 40}")

print np.dot(soaps_49492_1.get_soap(0,0), soaps_49492_1.get_soap(0,1999))
print "" * 20
print np.dot(soaps_49492_2.get_soap(0,0), soaps_49492_2.get_soap(0,1999))
print "" * 20
print np.dot(soaps_444412_1.get_soap(0,0), soaps_444412_1.get_soap(0,1999))
print "" * 20
print np.dot(soaps_49492_1.get_soap(0,1999), soaps_49492_2.get_soap(0,1999))
print "" * 20
print np.dot(soaps_49492_1.get_soap(0,1999), soaps_444412_1.get_soap(0,1999))
print "" * 20
print np.dot(soaps_49492_2.get_soap(0,1999), soaps_444412_1.get_soap(0,1999))

"""
soaps_49492_1 = soaplist.SOAPList("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr49Cu49Al2/c/md/lammps/100/1/extxyzf.xyz", verb = True)
#soaps_49492_1 = soaplist.SOAPList("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr49Cu49Al2/c/md/lammps/100/1/extxyzf_eos_1999.xyz", verb = True)
# on Al central atoms
soaps_49492_1.compute_average(cutoff = "3.75", 
                      l_max = "6", 
                      n_max = "8", 
                      n_Z = "1", 
                      Z = "{13}", 
                      n_species = "3", 
                      species_Z = "{13 29 40}")
#soaps_49492_1.write_avg_soaps()

nfms = soaps_49492_1.get_n_frames()
for ifms in range(nfms):
    print ifms, np.dot(soaps_49492_1.get_soap(0,nfms-1), soaps_49492_1.get_soap(0,ifms))

"""

