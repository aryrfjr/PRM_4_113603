#!/usr/bin/python
#

#
# This is a script that reads a <STEP>_<TRUE_STEP>_ICOHP.lobster 
# file generated with the script ML-ICOHP.py and creates a 
# per-interaction "data frame" as described on 07/08/2019-(5).
#

# Libraries
import numpy as np
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_cohp

INTERACTIONS = ["Al-Al", "Cu-Al", "Cu-Cu", "Zr-Al", "Zr-Cu", "Zr-Zr"]

WHAT_TO_DO = sys.argv[1]

if WHAT_TO_DO == "strains":
    ICOHP_DIR = sys.argv[2]
    STEPS_DUMP = ["1", "26", "51", "76", "101", "126", "151", "176", "201", "226", "251", "276", "401", "501"]
    STEPS_REAL = ["0", "500000", "1000000", "1500000", "2000000", "2500000", "3000000", "3500000", "4000000", "4500000", "5000000", "5500000", "8000000", "10000000"]
    STRAINS = ["0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0", "4.5", "5.0", "5.5", "8.0", "10.0"]
    ifiles = [None] * 6
    for ii in range(6):
        ifiles[ii] = open(ICOHP_DIR+"/"+INTERACTIONS[ii]+"_ICOHPLIST.csv","w")
        ifiles[ii].write("Strain,minusICOHP,A,B,Distance,TotInt,FracInt\n")

    for ifs in range(len(STEPS_DUMP)):
        info_file = ICOHP_DIR+"/"+STEPS_DUMP[ifs]+"_"+STEPS_REAL[ifs]+"_ICOHPLIST.info"
        fileobj = open(info_file)
        info_lines = fileobj.readlines()
        fileobj.close()
        tot_int = float(info_lines[0].split()[0]) # total number of interactions
        icohp_file = ICOHP_DIR+"/"+STEPS_DUMP[ifs]+"_"+STEPS_REAL[ifs]+"_ICOHPLIST.lobster"
        print icohp_file
        bonds = read_cohp(cohpfile = icohp_file, task = "ICOHPLIST")
        for ibk in range(len(bonds)):
            bk = bonds.keys()[ibk]
            b = bonds[bk]
            # Selecting the database
            sA = b.get_symbA()[:2]
            sB = b.get_symbB()[:2]
            try:
                ii = INTERACTIONS.index(sA+"-"+sB)
            except ValueError:
                ii = INTERACTIONS.index(sB+"-"+sA)
            count_int = float(info_lines[ii+3].split()[1])
            frac_int = 100.0 - (tot_int - count_int) / tot_int * 100.0
            ifiles[ii].write("\"%s\",%f,%s,%s,%f,%d,%f\n" % (STRAINS[ifs], -b.get_info(), 
                             b.get_symbA(), b.get_symbB(), b.get_distance(), int(count_int), frac_int))

elif WHAT_TO_DO == "ncs": # nominal compositions
    STEP_DUMP = "40"
    STEP_REAL = "101025000"
    STRAIN = "0.0"
    #X = [2, 6, 8, 10]
    #NCS = ["Zr51Cu39.5Al9.5", "Zr47Cu46Al7", "Zr40Cu54Al6"]
    NCS = ["Zr40Cu54Al6", "Zr47Cu46Al7", "Zr51Cu39.5Al9.5"]
    NCIDS = ["NC1", "NC2", "NC3"]
    WRITE_DIR = sys.argv[2]

    ifiles = [None] * 6
    for ii in range(6):
        ifiles[ii] = open(WRITE_DIR+"/"+INTERACTIONS[ii]+"_ICOHPLIST.csv","w")
        ifiles[ii].write("NC,minusICOHP,A,B,Distance,TotInt,FracInt\n")

    for inc in range(len(NCS)):
        #NC = "Zr"+str(50-(X[inc]/2))+"Cu"+str(50-(X[inc]/2))+"Al"+str(X[inc])
        NC = NCS[inc]
        #ICOHP_DIR = "/home/aryjr/UFSCar/MG-NMR/SS-ML/"+NC+"/c/md/lammps/ultimate/2/icohps_3"
        ICOHP_DIR = "/home/aryjr/UFSCar/MG-NMR/SS-ML/"+NC+"/icohps"
        info_file = ICOHP_DIR+"/"+STEP_DUMP+"_"+STEP_REAL+"_ICOHPLIST.info"
        fileobj = open(info_file)
        info_lines = fileobj.readlines()
        fileobj.close()
        tot_int = float(info_lines[0].split()[0]) # total number of interactions
        icohp_file = ICOHP_DIR+"/"+STEP_DUMP+"_"+STEP_REAL+"_ICOHPLIST.lobster"
        print icohp_file
        bonds = read_cohp(cohpfile = icohp_file, task = "ICOHPLIST")
        for ibk in range(len(bonds)):
            bk = bonds.keys()[ibk]
            b = bonds[bk]
            # Selecting the database
            sA = b.get_symbA()[:2]
            sB = b.get_symbB()[:2]
            try:
                ii = INTERACTIONS.index(sA+"-"+sB)
            except ValueError:
                ii = INTERACTIONS.index(sB+"-"+sA)
            count_int = float(info_lines[ii+3].split()[1])
            frac_int = 100.0 - (tot_int - count_int) / tot_int * 100.0
            ifiles[ii].write("\"%s\",%f,%s,%s,%f,%d,%f\n" % (NCIDS[inc], -b.get_info(), 
                             b.get_symbA(), b.get_symbB(), b.get_distance(), int(count_int), frac_int))

for ii in range(6):
    ifiles[ii].close()

