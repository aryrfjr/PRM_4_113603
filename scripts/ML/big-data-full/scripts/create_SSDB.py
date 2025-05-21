#!/usr/bin/python

#
# Creates a single-SOAP (SS) database for all <ID_RUN>s in a given 
# nominal composition. The code follows the directory structure 
# below (see on 05/10/2018-(10)).
# 
# <CHEM_COMPOSITION>/c/md/lammps/100/<ID_RUN>/2000/<SUB_STEP>
#
# with <SUB_STEP>s ranging from 0 to 14.
# 
# usage:
# 
# python create_SSDB.py <TYPE> <CHEM_COMPOSITION>
#
# with <TYPE> can be PC (per-cluster) or PB (per-bond).
#

# Libraries
import os
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_cohp
from ase.io import read
import numpy as np
import pickle as pkl

#
# Subroutines
#

def check_QE_ICOHP():
    # Proceeding with further checks
    statuses = []
    gen_statuses = []
    for ss in range(15):
        qe_file_exists = True
        qe_scf_ok = False
        qe_finished = False
        lob_file_exists = False
        icohp_bonds = -1
        lob_finished = False
        icohp_file_exists = False
        icohp_ok = False
        icohp_lines = -1
        general_status = DONE
        # Checking the QE calculation
        try:
            fileobj = open(RUN_DIR+"/2000/"+str(ss)+"/"+CHEM_COMPOSITION+".scf.out")
        except:
            qe_file_exists = False
        if qe_file_exists:
            lines = fileobj.readlines()
            fileobj.close()
            for i, line in enumerate(lines):
                if "convergence has been achieved in" in line: qe_scf_ok = True
                if "This run was terminated on" in line:
                    qe_finished = True
                    lob_file_exists = True # assuming that
        # Checking the LOBSTER calculation
        try:
            fileobj = open(RUN_DIR+"/2000/"+str(ss)+"/"+CHEM_COMPOSITION+".lb.out")
        except:
            lob_file_exists = False
        if lob_file_exists:
            icohp_file_exists = True # assuming that
            lines = fileobj.readlines()
            fileobj.close()
            for i, line in enumerate(lines):
                if "setting up CO interactions... found" in line:
                    icohp_bonds = int(line.split()[5])
                if "finished in" in line:
                    lob_finished = True
        # Checking the ICOHP results
        try:
            fileobj = open(RUN_DIR+"/2000/"+str(ss)+"/ICOHPLIST.lobster")
        except:
            icohp_file_exists = False
        if icohp_file_exists:
            lines = fileobj.readlines()
            fileobj.close()
            icohp_lines = len(lines) - 1
            icohp_ok = icohp_lines == icohp_bonds
        # setting the status
        if qe_file_exists and not qe_scf_ok:
            general_status = SCF_NOT_CONVERGED
        elif not lob_finished and icohp_lines > -1:
            general_status = LOB_NOT_FINISHED_BUT_ICOHP_OK
        elif qe_file_exists and qe_scf_ok and not lob_finished and icohp_lines > -1:
            general_status = SCF_CONVERGED_LOB_NOT_FINISHED_BUT_ICOHP_OK
        elif qe_file_exists and qe_scf_ok and icohp_lines == -1:
            general_status = SCF_CONVERGED_BUT_NO_ICOHP
        elif not qe_scf_ok and icohp_lines == -1:
            general_status = SCF_NOT_CONVERGED_AND_NO_ICOHP
        gen_statuses.append(general_status)
        statuses.append([qe_file_exists, qe_scf_ok, qe_finished, lob_file_exists, 
                         icohp_bonds, lob_finished, icohp_file_exists, icohp_ok, icohp_lines])
    return gen_statuses, statuses

def eval_cluster(bonds, acs, offset, cell, pos, symbs, bonds_to_writeK):
    symbA = symbs[acs[0]-1]
    knf = "KNF: " # keys not found
    sum_ICOHP = 0.0 # used only with a per-cluster database
    bonds_to_write = [] # used only with a per-bond database
    for ib in range(1, len(acs)):
        symbB = symbs[acs[ib]-1]
        key = symbA+str(acs[0])+"-"+symbB+str(acs[ib])
        nknf = ""
        try: # trying with the normal key AtomA-AtomB
            bond = bonds[key]
        except KeyError:
            nknf = key
            try: # trying with the inverted key AtomB-AtomA
                key = symbB+str(acs[ib])+"-"+symbA+str(acs[0])
                bond = bonds[key]
                nknf = ""
            except KeyError: # none of the two keys have been found
                nknf = key
        if nknf == "": # if at least one of the two keys have been found
            if TYPE == "PC": # a per-cluster database
                sum_ICOHP += -bond.get_info()
            elif TYPE == "PB": # a per-bond database
                if not bond.get_symbA()+"-"+bond.get_symbB() in bonds_to_writeK and not bond.get_symbB()+"-"+bond.get_symbA() in bonds_to_writeK:
                    bonds_to_writeK.append(bond.get_symbA()+"-"+bond.get_symbB())
                    bonds_to_write.append(bond)
        else:
            knf = knf + nknf + ", "
    if knf == "KNF: ": # ok! The LOBSTER and quippy clusters are compatible ...
        if TYPE == "PC": # a per-cluster database
            if symbA == "Zr":
                ftw = FDB_DONE_Zr # file to write
                stl = ALL_SOAPS_Zr # SOAPs to load
            elif symbA == "Cu":
                ftw = FDB_DONE_Cu
                stl = ALL_SOAPS_Cu
            elif symbA == "Al":
                ftw = FDB_DONE_Al
                stl = ALL_SOAPS_Al
            # ... write the .xyz file and ... 
            ftw.write("%d\n" % len(acs))
            ftw.write("Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\" ID_RUN=%d SUB_STEP=%d ID_CENTRAL_ATOM=%d SUM_ICOHP=%f\n" 
                      % (cell[0,0], cell[0,1], cell[0,2], cell[1,0], cell[1,1], cell[1,2], cell[2,0], cell[2,1], cell[2,2], int(id_run), sub_step, acs[0], sum_ICOHP))
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
            # ... saving the SOAP vectors of the central atom
            # naturally synchronized with the .xyz files written above
            stl.append(SS_SOAPS[acs[0]-1])
        elif TYPE == "PB": # a per-bond database
            for ibnd in range(len(bonds_to_write)):
                bond = bonds_to_write[ibnd]
                symbA = bond.get_symbA()[:2]
                symbB = bond.get_symbB()[:2]
                if symbA+symbB == "ZrZr":
                    ftw = FDB_DONE_ZrZr # file to write
                    stl = ALL_SOAPS_ZrZr # SOAPs to load
                elif symbA+symbB == "CuCu":
                    ftw = FDB_DONE_CuCu
                    stl = ALL_SOAPS_CuCu
                elif symbA+symbB == "AlAl":
                    ftw = FDB_DONE_AlAl
                    stl = ALL_SOAPS_AlAl
                elif symbA+symbB == "ZrCu" or symbA+symbB == "CuZr":
                    ftw = FDB_DONE_ZrCu
                    stl = ALL_SOAPS_ZrCu
                elif symbA+symbB == "ZrAl" or symbA+symbB == "AlZr":
                    ftw = FDB_DONE_ZrAl
                    stl = ALL_SOAPS_ZrAl
                elif symbA+symbB == "CuAl" or symbA+symbB == "AlCu":
                    ftw = FDB_DONE_CuAl
                    stl = ALL_SOAPS_CuAl
                # ... write the .bnd file and ... 
                ftw.write("%d %d %s %s %f %f\n" % (int(id_run), sub_step, 
                          int(bond.get_symbA()[2:]), int(bond.get_symbB()[2:]), bond.get_distance(), -bond.get_info()))
                # ... saving the SOAP vectors in a dictionary
                # to avoid repeated elements
                keysoap = id_run+"-"+str(sub_step)+"-"+bond.get_symbA()[2:]
                if not keysoap in stl.keys():
                    stl[keysoap] = SS_SOAPS[int(bond.get_symbA()[2:])-1]
                    #if stl == ALL_SOAPS_AlAl: print keysoap
                keysoap = id_run+"-"+str(sub_step)+"-"+bond.get_symbB()[2:]
                if not keysoap in stl.keys():
                    stl[keysoap] = SS_SOAPS[int(bond.get_symbB()[2:])-1]
                    #if stl == ALL_SOAPS_AlAl: print keysoap
    else:
        # ... otherwise, I will write the unmatched bonds (UBs)
        # and this cluster will be discarded even with a DONE status.
        # It is IMPORTANT to point that these are bonds that were 
        # detected only by quippy (see more on 20/06/2019-(6)). Bonds 
        # that were detected only by LOBSTER are naturally discarded 
        # since the quippy clusters are those with which the SOAP 
        # vectors are generated and I'm supposed to have the ICOHP 
        # values for all bonds in the clusters.
        FDB_DUB.write("%d %d %s\n" % (int(id_run), sub_step, knf[:len(knf)-2]))
    return bonds_to_writeK # used only with a per-bond database

def write_files():
    ssdir = RUN_DIR+"/2000/"+str(sub_step)+"/"
    # Reading ICOHPLIST.lobster
    bonds = read_cohp(cohpfile = ssdir+"ICOHPLIST.lobster", task = "ICOHPLIST")
    # Reading the .xyz file
    atoms = read(ssdir+CHEM_COMPOSITION+".xyz")
    cell = atoms.get_cell()
    pos = atoms.get_positions()
    symbs = atoms.get_chemical_symbols()
    # Reading the reference file with the clusters
    fileobj = open(ssdir+"lobsterin-quippy")
    lines = fileobj.readlines()
    fileobj.close()
    lastA = -1
    # Setting up each one
    bonds_to_writeK = [] # used only with a per-bond database
    for i in range(9, len(lines)):
        refA = int(lines[i].split()[2])
        refB = int(lines[i].split()[4])
        if refA != lastA:
            if i > 9: bonds_to_writeK = eval_cluster(bonds, acs, offset, cell, pos, symbs, bonds_to_writeK)
            acs = []
            acs.append(refA)
            offset = []
        acs.append(refB)
        offset.append([int(lines[i].split()[6]), int(lines[i].split()[7]), int(lines[i].split()[8])])
        lastA = refA
    # Setting up the last
    bonds_to_writeK = eval_cluster(bonds, acs, offset, cell, pos, symbs, bonds_to_writeK)

#
# The script
#

# Constants and arguments
BD_DIR = "/home/aryjr/UFSCar/MG-NMR/ML/big-data-full"
TYPE = sys.argv[1]
if TYPE != "PC" and TYPE != "PB":
    sys.exit("The <TYPE> argument must be PC (per-cluster) or PB (per-bond)!")
CHEM_COMPOSITION = sys.argv[2]
DB_DIR = BD_DIR+"/"+CHEM_COMPOSITION+"-"+TYPE+"SSDB"
RUN_DIR_NOT_PREPARED = 0
DONE = 1
SCF_NOT_CONVERGED = 2
LOB_NOT_FINISHED_BUT_ICOHP_OK = 3
SCF_CONVERGED_LOB_NOT_FINISHED_BUT_ICOHP_OK = 4
SCF_CONVERGED_BUT_NO_ICOHP = 5
SCF_NOT_CONVERGED_AND_NO_ICOHP = 6
DONE_BUT_UNMATCHED_BONDS = 7

if os.path.exists(DB_DIR):
    if TYPE == "PC":
        sys.exit("The database directory <CHEM_COMPOSITION>-PCSSDB already exists, please rename it!")
    elif TYPE == "PB":
        sys.exit("The database directory <CHEM_COMPOSITION>-PBSSDB already exists, please rename it!")
else:
    os.mkdir(DB_DIR)

# An information file with those environments in runs with DONE status
# but containig unmatched bonds
FDB_DUB = open(DB_DIR+"/DUB.info","w")
# An information file with those <SUB_STEP>s with a non-DONE status
FDB_ND = open(DB_DIR+"/ND.info","w")
if TYPE == "PC": # a per-cluster database
    # Per-element big "Extended xyz" files with all environments
    FDB_DONE_Zr = open(DB_DIR+"/Zr.xyz","w")
    FDB_DONE_Cu = open(DB_DIR+"/Cu.xyz","w")
    FDB_DONE_Al = open(DB_DIR+"/Al.xyz","w")
    # And also all SOAP vectors
    ALL_SOAPS_Zr = []
    ALL_SOAPS_Cu = []
    ALL_SOAPS_Al = []
elif TYPE == "PB": # a per-bond database
    # Per-bond-type .bnd big files with ICOHP information
    FDB_DONE_ZrZr = open(DB_DIR+"/Zr-Zr.bnd","w")
    FDB_DONE_CuCu = open(DB_DIR+"/Cu-Cu.bnd","w")
    FDB_DONE_AlAl = open(DB_DIR+"/Al-Al.bnd","w")
    FDB_DONE_ZrCu = open(DB_DIR+"/Zr-Cu.bnd","w")
    FDB_DONE_ZrAl = open(DB_DIR+"/Zr-Al.bnd","w")
    FDB_DONE_CuAl = open(DB_DIR+"/Cu-Al.bnd","w")
    # And also all SOAP vectors as dictionaries
    ALL_SOAPS_ZrZr = {}
    ALL_SOAPS_CuCu = {}
    ALL_SOAPS_AlAl = {}
    ALL_SOAPS_ZrCu = {}
    ALL_SOAPS_ZrAl = {}
    ALL_SOAPS_CuAl = {}
# This is a general portion of the script
SS_SOAPS = None # <SUB_STEP>s SOAPs
# Loop over all LAMMPS <RUN>s
for id_run in os.listdir(BD_DIR+"/"+CHEM_COMPOSITION+"/c/md/lammps/100"):
    RUN_DIR = BD_DIR+"/"+CHEM_COMPOSITION+"/c/md/lammps/100/"+id_run
    SOAPS_DIR = BD_DIR+"/"+CHEM_COMPOSITION+"-SOAPS/c/md/lammps/100/"+id_run
    if os.path.isdir(RUN_DIR):
        # checking its <SUB_STEP>s's statuses
        gsts, sts = check_QE_ICOHP()
        # now writing the files
        if TYPE == "PC": # a per-cluster database
            print "Writing .xyz files for <ID_RUN> = %d ..." % int(id_run)
        elif TYPE == "PB": # a per-bond database
            print "Writing .bnd files for <ID_RUN> = %d ..." % int(id_run)
        for sub_step in range(0, 15):
            if gsts[sub_step] == DONE or gsts[sub_step] == SCF_CONVERGED_LOB_NOT_FINISHED_BUT_ICOHP_OK:
                ftl = open(SOAPS_DIR+"/2000/"+str(sub_step)+"/SOAPS.vec","rb")
                SS_SOAPS = pkl.load(ftl)
                ftl.close()
                write_files()
            else: # the <SUB_STEP> has a non-DONE status
                #print "Ops!!! <ID_RUN> = %d, <SUB_STEP> = %d" % (int(id_run), sub_step)
                FDB_ND.write("%d %d %d\n" % (int(id_run), sub_step, gsts[sub_step]))
if TYPE == "PC": # a per-cluster database
    # ok, done with the big .xyz files
    FDB_DONE_Zr.close()
    FDB_DONE_Cu.close()
    FDB_DONE_Al.close()
    # and finally dumping the per-element SOAP vectors synchronized with the respective .xyz files
    print "Dumping the soap vectors ..."
    # Zr
    ftd = open(DB_DIR+"/Zr-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_Zr, ftd)
    ftd.close()
    # Cu
    ftd = open(DB_DIR+"/Cu-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_Cu, ftd)
    ftd.close()
    # Al
    ftd = open(DB_DIR+"/Al-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_Al, ftd)
    ftd.close()
elif TYPE == "PB": # a per-bond database
    # ok, done with the big .bnd files
    FDB_DONE_ZrZr.close()
    FDB_DONE_CuCu.close()
    FDB_DONE_AlAl.close()
    FDB_DONE_ZrCu.close()
    FDB_DONE_ZrAl.close()
    FDB_DONE_CuAl.close()
    # and finally dumping the per-element SOAP vectors synchronized with the respective .xyz files
    print "Dumping the soap vectors ..."
    ftd = open(DB_DIR+"/Zr-Zr-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_ZrZr, ftd)
    ftd.close()
    ftd = open(DB_DIR+"/Cu-Cu-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_CuCu, ftd)
    ftd.close()
    ftd = open(DB_DIR+"/Al-Al-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_AlAl, ftd)
    ftd.close()
    ftd = open(DB_DIR+"/Zr-Cu-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_ZrCu, ftd)
    ftd.close()
    ftd = open(DB_DIR+"/Zr-Al-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_ZrAl, ftd)
    ftd.close()
    ftd = open(DB_DIR+"/Cu-Al-SOAPS.vec","wb")
    pkl.dump(ALL_SOAPS_CuAl, ftd)
    ftd.close()
# Finishing
FDB_DUB.close()
FDB_ND.close()
print "Done!!!"

