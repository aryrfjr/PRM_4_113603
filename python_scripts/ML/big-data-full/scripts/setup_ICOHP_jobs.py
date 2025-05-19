#!/usr/bin/python

#
# This script will setup the Slurm job files 
# to run in parallel 30 ICOHP calculations, 
# 15 for each <ID_RUN>. Recalling that only 
# the last <STEP> will be considered, for 
# which 14 <SUB_STEP> will  be derived to 
# generated the distorted structures following 
# the following directory structure (see on 
# 05/10/2018-(10))
# 
# <CHEM_COMPOSITION>/<ID_RUN>/<STEP>/<SUB_STEP>
#
# usage:
# 
# python setup_ICOHP_jobs.py <CHEM_COMPOSITION> <ID_RUN>
#

# Libraries
import os
import quippy
# from quippy import Atoms # never used in this script and not found in current quippy
from quippy import descriptors
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.atom.desc.soap import soaplist
from theo4m.io import read_lammps
import numpy as np
import pickle as pkl

#
# Subroutines
#

# writing the extended .xyz file to be used by quippy
def dump_soaps(f2d, soaps):
    ftd = open(f2d,"wb")
    pkl.dump(soaps, ftd)
    ftd.close()

# writing the extended .xyz file to be used by quippy
def write_ext_xyz(f2w, tcell, atoms, spc_symbs):
    ftw = open(f2w,"w")
    ftw.write("%d\n" % len(atoms))
    # The order of the lattice vector was wrong in a previous version 
    # but not now it is correct. See more at:
    # https://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz
    ftw.write("Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" 
              % (tcell[0,0], tcell[0,1], tcell[0,2], tcell[1,0], tcell[1,1], tcell[1,2], tcell[2,0], tcell[2,1], tcell[2,2]))
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()
    geom = np.dot(geom, tcell)
    for i in range(len(atoms)):
        ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
    ftw.close()

# writing the .scf.in file to be used by QE
def write_qe_inp(f2w, tcell, atoms, prefix, outdir, spc_symbs, spc_nbfs, check = False):
    nbnd = 0
    symbs = atoms.get_chemical_symbols()
    for i in range(len(spc_symbs)):
        nbnd += symbs.count(spc_symbs[i]) * spc_nbfs[i]
    nbnd += int(nbnd * 0.1)
    ftw = open(f2w,"w")
    ftw.write(" &CONTROL\n")
    ftw.write("                       title = \'"+prefix+"\',\n")
    ftw.write("                 calculation = \'scf\',\n")
    ftw.write("                restart_mode = \'from_scratch\',\n")
    ftw.write("                      outdir = \'"+outdir+"\',\n")
    ftw.write("                      wfcdir = \'"+outdir+"\',\n")
    ftw.write("                  pseudo_dir = \'/scratch/rmnvm/ary.junior/UPFs/\',\n")
    ftw.write("                      prefix = \'"+prefix+"\',\n")
    ftw.write("                   verbosity = \'high\',\n")
    ftw.write("                     tprnfor = .true.,\n")
    ftw.write("                  wf_collect = .true.\n")
    ftw.write(" /\n")
    ftw.write(" &SYSTEM\n")
    ftw.write("                       ibrav = 0,\n")
    if check: ftw.write("                   celldm(1) = 1.8897265,\n")
    ftw.write("                         nat = %d,\n" % len(atoms))
    ftw.write("                        ntyp = 3,\n")
    ftw.write("                        nbnd = %d,\n" % nbnd)
    ftw.write("                     ecutwfc = 70,\n")
    ftw.write("                     ecutrho = 560,\n")
    ftw.write("                 occupations = \'smearing\',\n")
    ftw.write("                    smearing = \'fd\',\n")
    ftw.write("                     degauss = 0.008,\n")
    ftw.write("                       nosym = .true.,\n")
    ftw.write("                       noinv = .true.\n")
    ftw.write(" /\n")
    ftw.write(" &ELECTRONS\n")
    ftw.write("            electron_maxstep = 100,\n")
    ftw.write("                    conv_thr = 1.0D-6,\n")
    ftw.write("                 mixing_beta = 0.4\n")
    ftw.write(" /\n")
    if check:
        ftw.write("CELL_PARAMETERS alat\n")
    else:
        ftw.write("CELL_PARAMETERS {angstrom}\n")
    ftw.write("%10f %10f %10f\n" % (tcell[0,0], tcell[0,1], tcell[0,2]))
    ftw.write("%10f %10f %10f\n" % (tcell[1,0], tcell[1,1], tcell[1,2]))
    ftw.write("%10f %10f %10f\n" % (tcell[2,0], tcell[2,1], tcell[2,2]))
    ftw.write("ATOMIC_SPECIES\n")
    ftw.write("   Al  26.9815   Al.aryjr-1.0.0-ld1.UPF\n")
    ftw.write("   Cu  63.5463   Cu.aryjr-1.0.0-ld1.UPF\n")
    ftw.write("   Zr  91.2240   Zr.aryjr-1.0.0-ld1.UPF\n")
    ftw.write("ATOMIC_POSITIONS {crystal}\n")
    geom = atoms.get_positions()
    for i in range(len(atoms)):
        ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
    ftw.write("K_POINTS automatic\n")
    ftw.write("  2 2 2   1 1 1\n")
    ftw.close()

# Using quippy to create the lobsterin file
def write_lobsterin(f2w, cats, qatoms, auto=False, lower=0.0, higher=0.0):
    ftw = open(f2w,"w")
    ftw.write("skipDOS\n")
    ftw.write("skipCOOP\n")
    ftw.write("skipPopulationAnalysis\n")
    ftw.write("COHPstartEnergy -100.0\n")
    ftw.write("COHPendEnergy 40.0\n")
    ftw.write("basisSet Bunge\n")
    ftw.write("basisfunctions Al 3s 3p\n")
    ftw.write("basisfunctions Cu 3s 4s 3p 3d\n")
    ftw.write("basisfunctions Zr 4s 5s 4p 4d\n")
    if auto:
        ftw.write("cohpGenerator from %1.2f to %1.2f\n" % (lower, higher))
    else:
        csymbs = qatoms.get_chemical_symbols()
        # loop over central atoms
        for k in range(len(cats)):
            indices, offsets = qatoms.connect.get_neighbours(cats[k][0] + 1)
            for i, offset in zip(indices, offsets):
                ftw.write("cohpbetween atom %d atom %d offset %d %d %d\n" % (cats[k][0] + 1, i, offset[0], offset[1], offset[2]))
    ftw.close()

# Shears a cell in a single shear matrix parameter
def shear_cell(a, b, c, par, value):
    # http://www.cs.brandeis.edu/~cs155/Lecture_07_6.pdf
    # https://www.mauriciopoppe.com/notes/computer-graphics/transformation-matrices/shearing/
    ash = bsh = csh = dsh = esh = fsh = 0.0
    # Parameters
    # 0 ==> ash
    # 1 ==> bsh
    # 2 ==> csh
    # 3 ==> dsh
    # 4 ==> esh
    # 5 ==> fsh
    if par == 0:
        ash = value
    elif par == 1:
        bsh = value
    elif par == 2:
        csh = value
    elif par == 3:
        dsh = value
    elif par == 4:
        esh = value
    elif par == 5:
        fsh = value
    # The shearing array
    shear = np.array([[1.0, ash, bsh, 0.0],
                      [csh, 1.0, dsh, 0.0],
                      [esh, fsh, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0]])
    # The 
    lattv = np.array([[a, 0.0, 0.0],
                      [0.0, b, 0.0],
                      [0.0, 0.0, c],
                      [1.0, 1.0, 1.0]])
    tcell = np.matmul(shear, lattv)
    return np.array([[tcell[0,0], tcell[0,1], tcell[0,2]],
                     [tcell[1,0], tcell[1,1], tcell[1,2]],
                     [tcell[2,0], tcell[2,1], tcell[2,2]]])

# Compress or tension a cell along x, y, or z
def comp_tens_cell(a, b, c, par, value_comp, value_tens):
    # http://www.cs.brandeis.edu/~cs155/Lecture_07_6.pdf
    # https://www.mauriciopoppe.com/notes/computer-graphics/transformation-matrices/shearing/
    ash = bsh = csh = dsh = esh = fsh = 0.0
    # Parameters
    # 0 ==> tens x
    # 1 ==> comp x
    # 2 ==> tens y
    # 3 ==> comp y
    # 4 ==> tens z
    # 5 ==> comp z
    # 6 ==> tens x, y, and z
    # 7 ==> comp x, y, and z
    if par == 0:
        a += a * value_tens
    elif par == 1:
        a -= a * value_comp
    elif par == 2:
        b += b * value_tens
    elif par == 3:
        b -= b * value_comp
    elif par == 4:
        c += c * value_tens
    elif par == 5:
        c -= c * value_comp
    elif par == 6:
        a += a * value_tens
        b += b * value_tens
        c += c * value_tens
    elif par == 7:
        a -= a * value_comp
        b -= b * value_comp
        c -= c * value_comp
    return np.array([[a, 0.0, 0.0],
                      [0.0, b, 0.0],
                      [0.0, 0.0, c]])

#
# The script
#

# Constants
BD_DIR = "PRM_4_113603_repr/ML/big-data-full/"
LOCAL_DIR = "/home/aryjr/MLOps/"
SD_SCRATCH = "/scratch/rmnvm/ary.junior/"
REF_SYMBS = ["Zr","Cu","Al"] # the output must be in the order of LAMMPS
Z_NBFS = [10, 10, 4] # number of basis functions in the .UPF and in lobsterin to compute nbnd

# Input arguments
nominal_comp = sys.argv[1] # e.g. Zr49Cu49Al2
id_run = sys.argv[2] # e.g. 21
# Loop over the runs
#
# The first step is to load the reference 
# cell for the current <ID_RUN_1> (the last 
# MD step 2000) as a Theo4M AtomicSites object 
# so that I can generate the other 11 sub-steps
#
# Note: 2000 is the fixed <STEP> but in some 
# corrupted files (see on 20/06/2019-(3)) it 
# may be not and therefore I'm taking the last 
# one (-2).
id_run_dir = LOCAL_DIR+BD_DIR+nominal_comp+"/c/md/lammps/100/"+id_run
rsteps, ratoms, rccell = read_lammps(lmpoutput = id_run_dir+"/zca-th300.dump", 
                                 spc_symbs = REF_SYMBS, 
                                 frac = False, items = [-2])
atoms = ratoms[0]
ccell = rccell[0]
a = ccell[0,0]
b = ccell[1,1]
c = ccell[2,2]
# The original cell
# <SUB_STEP>: 0
sub_step = 0
outdir = SD_SCRATCH+"work/"+BD_DIR+nominal_comp+"/c/md/lammps/100/"+id_run+"/2000"
jobdir = LOCAL_DIR+BD_DIR+nominal_comp+"/c/md/lammps/100/"+id_run+"/2000"
soapsdir = LOCAL_DIR+BD_DIR+nominal_comp+"-SOAPS/c/md/lammps/100/"+id_run+"/2000"
if not os.path.exists(jobdir): os.mkdir(jobdir)
if not os.path.exists(soapsdir): os.makedirs(soapsdir)
# Writing the .scf.in file to be used by QE
if not os.path.exists(jobdir+"/"+str(sub_step)): os.mkdir(jobdir+"/"+str(sub_step))
if not os.path.exists(soapsdir+"/"+str(sub_step)): os.mkdir(soapsdir+"/"+str(sub_step))
write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+"-check.scf.in", ccell, 
             atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS, check = True)
write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".scf.in", ccell, 
             atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS)
write_ext_xyz(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".xyz", ccell, atoms, REF_SYMBS)
sub_step += 1
# Now the transformations (see on 04/01/2019-(3)):
# 1) shearing (angle to be defined) along the x, y, and z axis (6 new structures); 
for si in range(6): # <SUB_STEP>: 1 2 3 4 5 6
    # Applying the transformation to the cell
    scell = shear_cell(a, b, c, si, 0.5)
    # Writing the .scf.in file to be used by QE
    if not os.path.exists(jobdir+"/"+str(sub_step)): os.mkdir(jobdir+"/"+str(sub_step))
    if not os.path.exists(soapsdir+"/"+str(sub_step)): os.mkdir(soapsdir+"/"+str(sub_step))
    write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+"-check.scf.in", scell, 
                 atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS, check = True)
    write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".scf.in", scell, 
                 atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS)
    write_ext_xyz(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".xyz", scell, atoms, REF_SYMBS)
    sub_step += 1
# 2) ?% of compression and ?% of tension along the x, y, and z axis (6 new structures); 
for si in range(6): # <SUB_STEP>: 7 8 9 10 11 12
    ctcell = comp_tens_cell(a, b, c, si, 0.2, 0.02)
    # Writing the .scf.in file to be used by QE
    if not os.path.exists(jobdir+"/"+str(sub_step)): os.mkdir(jobdir+"/"+str(sub_step))
    if not os.path.exists(soapsdir+"/"+str(sub_step)): os.mkdir(soapsdir+"/"+str(sub_step))
    write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+"-check.scf.in", ctcell, 
                 atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS, check = True)
    write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".scf.in", ctcell, 
                 atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS)
    write_ext_xyz(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".xyz", ctcell, atoms, REF_SYMBS)
    sub_step += 1
# 3) ?% of isotropic compression and ?% of tension (2 new structures)
for si in range(6, 8): # <SUB_STEP>: 13 14
    ictcell = comp_tens_cell(a, b, c, si, 0.1, 0.01)
    # Writing the .scf.in file to be used by QE
    if not os.path.exists(jobdir+"/"+str(sub_step)): os.mkdir(jobdir+"/"+str(sub_step))
    if not os.path.exists(soapsdir+"/"+str(sub_step)): os.mkdir(soapsdir+"/"+str(sub_step))
    write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+"-check.scf.in", ictcell, 
                 atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS, check = True)
    write_qe_inp(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".scf.in", ictcell, 
                 atoms, nominal_comp, outdir+"/"+str(sub_step), REF_SYMBS, Z_NBFS)
    write_ext_xyz(jobdir+"/"+str(sub_step)+"/"+nominal_comp+".xyz", ictcell, atoms, REF_SYMBS)
    sub_step += 1
# Next, I can compute the SOAPS of the 15 cells
# and generate the files for the ICOHP calculations
# sys.setrecursionlimit(17000) # may be necessary
for iss in range(15):
    xyz_file = jobdir+"/"+str(iss)+"/"+nominal_comp+".xyz"
    # Firstly the single-SOAPs
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
    # Only the last step will be used
    soaps, cats, frame = soapsl.get_pasoaps(0,0)
    #if iss == 0:
    #    symbs = atoms.get_chemical_symbols()
    #    for i in range(0,len(soaps)):
    #        print symbs[cats[2][0]], symbs[cats[i][0]], np.dot(soaps[2], soaps[i])
    # dumping the soap vectors
    dump_soaps(soapsdir+"/"+str(iss)+"/SOAPS.vec", soaps)
    # Writing the lobsterin file
    write_lobsterin(jobdir+"/"+str(iss)+"/lobsterin-quippy", cats, frame)
    # distances based on the g(r) from my Sci. Rep.
    write_lobsterin(jobdir+"/"+str(iss)+"/lobsterin", cats, frame, auto=True, lower=1.8, higher=4.2)
