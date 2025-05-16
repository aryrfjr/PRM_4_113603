#!/usr/bin/python

import numpy as np
import quippy
import os, sys
sys.path.append("/home/aryjr/dev-python/my_lib")
from aseext.om.io import read_lammps
from quippy import Atoms
from quippy import descriptors

#
def create_desc(at, argsd):
    # only Al local environments
    desc = descriptors.Descriptor(argsd)
    at.set_cutoff(desc.cutoff())
    at.calc_connect()
    soap = desc.calc(at, grad=False)["descriptor"]
    # TODO: for some reason, average=False results in a weird vector here!
    centers = desc.calc(at, grad=False)["descriptor_index_0based"]
    return soap, centers

# writing the .xyz file to be used by quippy
def write_quippy_xyz(f2w, a, b, c, atoms):
    ftw = open(f2w,"w")
    ftw.write("%d\n" % len(atoms))
    ftw.write("Lattice=\"%10f 0.0 0.0 0.0 %10f 0.0 0.0 0.0 %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" % (a, b, c))
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()
    for i in range(len(atoms)):
        ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
    ftw.close()

# writing the .scf.in file to be used by QE
def write_qe_inp(f2w, a, b, c, atoms, prefix, outdir, spc_symbs, spc_nbfs):
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
    ftw.write("            electron_maxstep = 200,\n")
    ftw.write("                    conv_thr = 1.0D-6,\n")
    ftw.write("                 mixing_beta = 0.5\n")
    ftw.write(" /\n")
    ftw.write("CELL_PARAMETERS {angstrom}\n")
    ftw.write("%10f %10f %10f\n" % (a, 0.0, 0.0))
    ftw.write("%10f %10f %10f\n" % (0.0, b, 0.0))
    ftw.write("%10f %10f %10f\n" % (0.0, 0.0, c))
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

# using quippy to create the lobsterin file
def write_lobsterin(f2w, catis, qatoms):
    ftw = open(f2w,"w")
    ftw.write("COHPstartEnergy -100.0\n")
    ftw.write("COHPendEnergy 40.0\n")
    ftw.write("basisSet Bunge\n")
    ftw.write("basisfunctions Al 3s 3p\n")
    ftw.write("basisfunctions Cu 3s 4s 3p 3d\n")
    ftw.write("basisfunctions Zr 4s 5s 4p 4d\n")
    csymbs = qatoms.get_chemical_symbols()
    for k in range(len(catis)):
        indices, offsets = qatoms.connect.get_neighbours(catis[k] + 1)
        for i, offset in zip(indices, offsets):
            ftw.write("cohpbetween atom %d atom %d\n" % (catis[k] + 1, i))
    ftw.close()

# writing the slurm file to submit the QE scf calculation
def write_slurm_qe_scf(f2w, prefix, nseq, homesddir, basebddir, jobdir):
    ftw = open(f2w,"w")
    ftw.write("#!/bin/bash\n")
    ftw.write("#SBATCH --nodes=8                   # Numero de Nos\n")
    ftw.write("#SBATCH --ntasks-per-node=24        # Numero de tarefas por No\n")
    ftw.write("#SBATCH --ntasks=192                # Numero total de tarefas MPI\n")
    ftw.write("#SBATCH -p cpu_long                # Fila (partition) a ser utilizada\n")
    ftw.write("#SBATCH -J ZCA-qs-%d                # Nome job\n" % nseq)
    ftw.write("#SBATCH --exclusive                 # Utilizacao exclusiva dos nos durante a execucao do job\n")
    ftw.write("#SBATCH -o stdout.%j                # File to which STDOUT will be written\n")
    ftw.write("#SBATCH -e stderr.%j                # File to which STDERR will be written\n")
    ftw.write("#SBATCH --mail-type=END             # Type of email notification- BEGIN,END,FAIL,ALL\n")
    ftw.write("#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent\n")
    ftw.write("#SBATCH --time=05:00:00             # Altera o tempo limite HH:MM:SS\n")
    ftw.write("\n")
    ftw.write("# Exibe os nos alocados para o Job\n")
    ftw.write("echo $SLURM_JOB_NODELIST\n")
    ftw.write("nodeset -e $SLURM_JOB_NODELIST\n")
    ftw.write("\n")
    ftw.write("#\n")
    ftw.write("# MPI and Multithreads\n")
    ftw.write("#\n")
    ftw.write("export OMP_NUM_THREADS=1\n")
    ftw.write("export MP_TASK_AFFINITY=cpu:1\n")
    ftw.write("export MP_SHARED_MEMORY=yes\n")
    ftw.write("export MEMORY_AFFINITY=MCM\n")
    ftw.write("\n")
    ftw.write("#\n")
    ftw.write("# Modulos\n")
    ftw.write("#\n")
    ftw.write("source /scratch/app/modulos/intel-psxe-2017.sh\n")
    ftw.write("module load fftw/3.3.5_intel\n")
    ftw.write("export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so\n")
    ftw.write("\n")
    ftw.write("#\n")
    ftw.write("# Command lines\n")
    ftw.write("#\n")
    ftw.write("export CMD_PW=\"srun -n 192 $SCRATCH/opt/qe-6.2_ARY_FFTW3_EXT/bin/pw.x -npool 8\"\n")
    ftw.write("\n")
    ftw.write("export BASE_DIR=\""+basebddir+"/"+jobdir+"\"\n")
    ftw.write("export WORK_DIR=\"$SCRATCH/work/$BASE_DIR\"\n")
    ftw.write("export RUN_DIR=\"$SCRATCH/$BASE_DIR\"\n")
    ftw.write("\n")
    ftw.write("mkdir -p $WORK_DIR\n")
    ftw.write("\n")
    ftw.write("chmod -R 760 $RUN_DIR\n")
    ftw.write("chmod -R 760 $WORK_DIR\n")
    ftw.write("\n")
    ftw.write("cd $RUN_DIR\n")
    ftw.write("\n")
    ftw.write("$CMD_PW < $RUN_DIR/"+prefix+".scf.in > $RUN_DIR/"+prefix+".scf.out\n")
    ftw.write("\n")
    ftw.write("chmod -R 760 $RUN_DIR\n")
    ftw.close()

# writing the slurm file to submit the LOBSTER calculation
def write_slurm_lobster(f2w, prefix, nseq, homesddir, basebddir, jobdir):
    ftw = open(f2w,"w")
    ftw.write("#!/bin/bash\n")
    ftw.write("#SBATCH --nodes=1               # Numero de Nos\n")
    ftw.write("#SBATCH --ntasks-per-node=1     # Numero de tarefas por No\n")
    ftw.write("#SBATCH --ntasks=1              # Numero total de tarefas MPI\n")
    ftw.write("#SBATCH --cpus-per-task=24      # Numero de threads\n")
    ftw.write("#SBATCH -p cpu_small            # Fila (partition) a ser utilizada\n")
    ftw.write("#SBATCH -J ZCA-lc-%d            # Nome job\n" % nseq)
    ftw.write("#SBATCH --exclusive             # Utilizacao exclusiva dos nos durante a execucao do job\n")
    ftw.write("#SBATCH -o stdout.%j            # File to which STDOUT will be written\n")
    ftw.write("#SBATCH -e stderr.%j            # File to which STDERR will be written\n")
    ftw.write("#SBATCH --mail-type=END         # Type of email notification- BEGIN,END,FAIL,ALL\n")
    ftw.write("#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent\n")
    ftw.write("#SBATCH --time=02:00:00         # Altera o tempo limite HH:MM:SS\n")
    ftw.write("\n")
    ftw.write("# Exibe os nos alocados para o Job\n")
    ftw.write("echo $SLURM_JOB_NODELIST\n")
    ftw.write("nodeset -e $SLURM_JOB_NODELIST\n")
    ftw.write("\n")
    ftw.write("#\n")
    ftw.write("# MPI and Multithreads\n")
    ftw.write("#\n")
    ftw.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
    ftw.write("export MP_TASK_AFFINITY=cpu:1\n")
    ftw.write("export MP_SHARED_MEMORY=yes\n")
    ftw.write("export MEMORY_AFFINITY=MCM\n")
    ftw.write("\n")
    ftw.write("#\n")
    ftw.write("# Modulos\n")
    ftw.write("#\n")
    ftw.write("source /scratch/app/modulos/intel-psxe-2017.sh\n")
    ftw.write("module load fftw/3.3.5_intel\n")
    ftw.write("export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so\n")
    ftw.write("\n")
    ftw.write("#\n")
    ftw.write("# Command lines\n")
    ftw.write("#\n")
    ftw.write("export CMD_LB=\"srun -n 1 -c $SLURM_CPUS_PER_TASK $SCRATCH/opt/lobster-3.0.0/lobster-3.0.0\"\n")
    ftw.write("\n")
    ftw.write("export BASE_DIR=\""+basebddir+"/"+jobdir+"\"\n")
    ftw.write("export WORK_DIR=\"$SCRATCH/work/$BASE_DIR\"\n")
    ftw.write("export RUN_DIR=\"$SCRATCH/$BASE_DIR\"\n")
    ftw.write("\n")
    ftw.write("cd $RUN_DIR\n")
    ftw.write("\n")
    ftw.write("$CMD_LB > $RUN_DIR/"+prefix+".lb.out\n")
    ftw.write("\n")
    ftw.write("chmod -R 760 $RUN_DIR\n")
    ftw.close()

########################################################################
# the script
########################################################################
# which atoms (same order of the create_atoms command)
spc_symbs = ["Zr","Cu","Al"] # the output must be in the order of LAMMPS
spc_nbfs = [10, 10, 4] # number of basis functions in the .UPF and in lobsterin to compute nbnd

# some settings
homesddir = "/scratch/rmnvm/ary.junior"
basebddir = "MG-NMR/ML/misc/4th_test_quippy"
lmax = 4 # spherical harmonics basis band limit
nmax = 5 # number of radial basis functions

# reading a set of jobs to be submitted
fileobj = open(sys.argv[1])
lines = fileobj.readlines()
fileobj.close()
iitem = 1
for i, line in enumerate(lines):
    if i > 0: # skipping the head
        chemc = line.split()[0]
        idrun = line.split()[1]
        step = line.split()[2]
        sstep = line.split()[3]
        # creating the directory
        jobdir = chemc+"/"+idrun+"/"+step+"/"+sstep
        print "creating: \'"+jobdir+"\'" 
        os.makedirs(jobdir, 0755);
        # reading the LAMMPS output file
        lof = "/home/aryjr/UFSCar/MG-NMR/BMG/"+chemc+"/c/md/lammps/2/"+idrun
        atoms, cell, ccell = read_lammps("DUMP_ATOM", lof, int(step), spc_symbs, 0, False, False, None)
        a = ccell[0][0]
        b = ccell[1][1]
        c = ccell[2][2]
        # writing the .xyz file to be used by quippy
        write_quippy_xyz(jobdir+"/"+chemc+".xyz", a, b, c, atoms)
        # reading the LAMMPS output file again
        # TODO: I shouldn't read again to get the crystal coordinates
        atoms, cell, ccell = read_lammps("DUMP_ATOM", lof, int(step), spc_symbs, 0, True, False, None)
        a = ccell[0][0]
        b = ccell[1][1]
        c = ccell[2][2]
        # writing the .scf.in file to be used by QE
        outdir = homesddir+"/work/"+basebddir+"/"+jobdir
        write_qe_inp(jobdir+"/"+chemc+".scf.in", a, b, c, atoms, chemc, outdir, spc_symbs, spc_nbfs)
        # using quippy to create the lobsterin file
        # cutoff from RDFs reported in [PRL 102, 245501 (2009)]
        soap_args = "soap cutoff=3.7 l_max="+str(lmax)+" n_max="+str(nmax)+" atom_sigma=0.5 n_Z=3 Z={13 29 40} n_species=3 species_Z={13 29 40} average=True"
        atl = quippy.AtomsList(jobdir+"/"+chemc+".xyz")
        soap, cats = create_desc(atl[len(atl)-1], soap_args)
        qatoms = atl[len(atl)-1]
        write_lobsterin(jobdir+"/lobsterin", cats[0], qatoms)
        # writing the slurm file to submit the QE scf calculation
        write_slurm_qe_scf(jobdir+"/ZCA-qs.sh", chemc, i, homesddir, basebddir, jobdir)
        # writing the slurm file to submit the LOBSTER calculation
        #write_slurm_lobster(jobdir+"/ZCA-lc.sh", chemc, i, homesddir, basebddir, jobdir)

