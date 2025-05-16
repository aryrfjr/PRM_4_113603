#!/usr/bin/python

#
# Check which steps and sub-steps of a 
# given run have been finished according 
# to the following directory structure (see 
# on 05/10/2018-(10))
# 
# <CHEM_COMPOSITION>/c/md/lammps/100/<ID_RUN>/<STEP>/<SUB_STEP>
#
# usage:
# 
# python check_QE_ICOHP.py <CHEM_COMPOSITION> <ID_RUN>
#

# Libraries
import sys
import os.path

#
# Subroutines
#

#
# The script
#

# Constants
BD_DIR = "MG-NMR/ML/big-data-full/"
LOCAL_DIR = "/home/aryjr/UFSCar/"

# Input arguments
nominal_comp = sys.argv[1]
id_run = sys.argv[2]
step = 2000 # only the last one for now

# Check if the directory has been prepared
if not os.path.exists(LOCAL_DIR+BD_DIR+nominal_comp+"/c/md/lammps/100/"+str(id_run)):
    sys.exit("ID_RUN_DIR_DOES_NOT_EXIST")

# Check if the directory has been prepared
RUN_DIR = LOCAL_DIR+BD_DIR+nominal_comp+"/c/md/lammps/100/"+str(id_run)+"/"+str(step)
if not os.path.exists(RUN_DIR):
    sys.exit("RUN_DIR_NOT_PREPARED")

# Proceeding with further checks
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
    general_status = "DONE"
    # Checking the QE calculation
    try:
        fileobj = open(RUN_DIR+"/"+str(ss)+"/"+nominal_comp+".scf.out")
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
        fileobj = open(RUN_DIR+"/"+str(ss)+"/"+nominal_comp+".lb.out")
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
        fileobj = open(RUN_DIR+"/"+str(ss)+"/ICOHPLIST.lobster")
    except:
        icohp_file_exists = False
    if icohp_file_exists:
        lines = fileobj.readlines()
        fileobj.close()
        icohp_lines = len(lines) - 1
        icohp_ok = icohp_lines == icohp_bonds
    # setting the status
    if qe_file_exists and not qe_scf_ok:
        general_status = "SCF_NOT_CONVERGED"
    elif not lob_finished and icohp_lines > -1:
        general_status = "LOB_NOT_FINISHED_BUT_ICOHP_OK"
    elif qe_file_exists and qe_scf_ok and not lob_finished and icohp_lines > -1:
        general_status = "SCF_CONVERGED_LOB_NOT_FINISHED_BUT_ICOHP_OK"
    elif qe_file_exists and qe_scf_ok and icohp_lines == -1:
        general_status = "SCF_CONVERGED_BUT_NO_ICOHP"
    elif not qe_scf_ok and icohp_lines == -1:
        general_status = "SCF_NOT_CONVERGED_AND_NO_ICOHP"
    # output
    print id_run, ss, qe_file_exists, qe_scf_ok, qe_finished, lob_file_exists, icohp_bonds, lob_finished, icohp_file_exists, icohp_ok, icohp_lines, general_status

