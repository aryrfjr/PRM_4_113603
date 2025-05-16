#!/bin/bash
#SBATCH --nodes=21                   # Numero de Nos
#SBATCH --ntasks-per-node=24        # Numero de tarefas por No
#SBATCH --ntasks=504                # Numero total de tarefas MPI
#SBATCH -p cpu                # Fila (partition) a ser utilizada
#SBATCH -J zca-grp                # Nome job
#SBATCH --exclusive                 # Utilizacao exclusiva dos nos durante a execucao do job
#SBATCH -o stdout.%j                # File to which STDOUT will be written
#SBATCH -e stderr.%j                # File to which STDERR will be written
#SBATCH --mail-type=END             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent
#SBATCH --time=02:00:00             # Altera o tempo limite HH:MM:SS

#Exibe os nós alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

echo "-----------------------------------------"
echo "Inicio do job:" `date`

#
# MPI and Multithreads
#
export OMP_NUM_THREADS=1
export MP_TASK_AFFINITY=cpu:1
export MP_SHARED_MEMORY=yes
export MEMORY_AFFINITY=MCM

#
# Modulos
#
source /scratch/app/modulos/intel-psxe-2017.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/app/mathlibs/fftw/3.3.8_openmpi-2.0_gnu/lib
module load openmpi/gnu/2.0.4.2 # after SD update on Feb 2019

### Definir variaveis para processamento com MPI da SGI ###

export MPI_MEMMAP_OFF=1
export MPI_GROUP_MAX=1000
export MPI_COMM_MAX=1000
export MPI_TYPE_MAX=100000

### A opcao -np indica o numero total de processos.  ###
export LMP_CMD="srun -n 504 /scratch/rmnvm/ary.junior/opt/lammps-16Feb16_ARY/src/lmp_mpi_fftw3_SD_update"

### Variaveis de controle da simulacao ###
export ID_RUN="${1}"
export SYS_NAME="${2}"
export BASE_DIR="/scratch/rmnvm/ary.junior/MG-NMR/SS-ML"
export G_ID="$ID_RUN"

### Rodando a simulacao ###
export G_DIR="$BASE_DIR/$SYS_NAME/c/md/lammps/ultimate/2/${3}"
export RUN_DIR="$G_DIR/$ID_RUN"
mkdir -p $RUN_DIR

export GONE="`cat $G_DIR/1_0_G1_S$G_ID.grp`"
export GTWO="`cat $G_DIR/1_0_G2_S$G_ID.grp`"

cp $BASE_DIR/ZrCuAl-2011.lammps.eam $RUN_DIR
cp $BASE_DIR/$SYS_NAME/$SYS_NAME-10000-1-1-1.data $RUN_DIR
cd $RUN_DIR

cat>$RUN_DIR/$SYS_NAME.lmp.inp<<EOF
# General setup
units metal
boundary p p p
atom_style atomic

# Reading the input geometry
read_data $SYS_NAME-10000-1-1-1.data

# The EAM potential
pair_style eam/alloy
# For all pairs of atom types (* *), reads the EAM potential file 
# ZrCuAl.eam.alloy. See more at
# http://lammps.sandia.gov/doc/pair_eam.html
# Note: the masses and the cutoffs are already specified in the EAM 
#       potential file.
# Note: the EAM file above is in the DYNAMO multi-element setfl format.
# Note: the elements in the EAM file are mapped to LAMMPS atom types 
#       specified with the create_atoms command above by specifying N 
#       additional arguments after the filename.
pair_coeff * * ZrCuAl-2011.lammps.eam Zr Cu Al

###########################################################################
#
# Equilibration of the loaded atoms at 300 K
#
###########################################################################

# a time step of 1 fs
timestep 0.001

# setting the velocities of all atoms as an uniform distribution
velocity all create 300 12345 mom yes rot no dist gaussian

#
# Following the link below I start with the NVT ensemble
#
# https://lammps.sandia.gov/threads/msg08725.html
#

# setting time integration
#
fix 1 all nvt temp 300 300 1

# set thermo output
thermo 1000
thermo_style custom step lx ly lz press pxx pyy pzz pe temp

# run for 100 picoseconds
run 100000

write_restart zca-equ-nvt-300.restart

unfix 1

###########################################################################
#
# Equilibration of the loaded atoms in NVE
#
###########################################################################

#
# And then the NVE ensemble to see if the temperature 
# stays constant.  If so, then the system is at equilibrium.

# setting time integration
#
fix 1 all nve

# run for 150 picoseconds
run 150000

write_restart zca-equ-nve.restart

# store final cell length for strain calculations
variable tmp equal "lx"
variable L0 equal \${tmp}
print "Initial Length, L0: \${L0}"

unfix 1

###########################################################################
#
# Input file for uniaxial compressive loading
#
###########################################################################

reset_timestep 0

# Grouping
group g1 id $GONE
group g2 id $GTWO

fix 1 all npt temp 300 300 1 y 0 0 1 z 0 0 1 drag 1

# this is the strain rate in s^-1
variable srate equal 1.0e7
# pointing that metal unit if time is ps (10^-12 s)
variable srate1 equal "-v_srate / 1.0e12"
fix 2 all deform 1 x erate \${srate1} units box remap x

# output strain and stress info to file
# for units metal, pressure is in [bars] = 100 [kPa] = 1/10000 [GPa]
# p2, p3, p4 are in GPa
variable strain equal "(lx - v_L0)/v_L0"
variable p1 equal "v_strain"
variable p2 equal "-pxx/10000"
variable p3 equal "-pyy/10000"
variable p4 equal "-pzz/10000"
fix def1 all print 100000 "\${p1} \${p2} \${p3} \${p4}" file comp_100.def1.txt screen no

# display thermo
thermo 20000
thermo_style custom step v_strain temp v_p2 v_p3 v_p4 ke pe press vol lx ly lz

# dumpinng the structures with the respective stress tensors
compute spa all stress/atom NULL
dump dspa all custom 20000 zca-dspa.dump id type x y z c_spa[1] c_spa[2] c_spa[3] c_spa[4] c_spa[5] c_spa[6]

# Compute mean-squared displacements (MSD) for each type of atom
compute g1c g1 msd
compute g2c g2 msd

fix g1f g1 ave/time 1 1 1 c_g1c[1] c_g1c[2] c_g1c[3] c_g1c[4] file G1.msd
fix g2f g2 ave/time 1 1 1 c_g2c[1] c_g2c[2] c_g2c[3] c_g2c[4] file G2.msd

run 20000000
EOF

$LMP_CMD < $RUN_DIR/$SYS_NAME.lmp.inp > $RUN_DIR/$SYS_NAME.lmp.out

echo "Final do job:" `date`
echo "-----------------------------------------" 

