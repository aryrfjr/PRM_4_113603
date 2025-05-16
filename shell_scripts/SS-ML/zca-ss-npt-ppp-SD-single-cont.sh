#!/bin/bash
#SBATCH --nodes=8                   # Numero de Nos
#SBATCH --ntasks-per-node=24        # Numero de tarefas por No
#SBATCH --ntasks=192                # Numero total de tarefas MPI
#SBATCH -p cpu_small                # Fila (partition) a ser utilizada
#SBATCH -J zca-ss-c                # Nome job
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
export LMP_CMD="srun -n 192 /scratch/rmnvm/ary.junior/opt/lammps-16Feb16_ARY/src/lmp_mpi_fftw3_SD_update"
export BASE_DIR="/scratch/rmnvm/ary.junior/MG-NMR/SS-ML"

### Variaveis de controle da simulacao ###
export SYS_NAME="${1}"
export ID_RUN_RESTART=${2}

export ID_RUN="${3}"
export SRATE="${4}"
export N_THERMO="${5}"
export N_RUN="${6}"

### Rodando a simulacao ###
export RUN_DIR="$BASE_DIR/$SYS_NAME/c/md/lammps/$ID_RUN-"
mkdir -p $RUN_DIR

cp $BASE_DIR/ZrCuAl-2011.lammps.eam $RUN_DIR
cp $BASE_DIR/$SYS_NAME/c/md/lammps/$ID_RUN_RESTART/zca-ss-npt.restart $RUN_DIR/zca-ss-npt_PREV.restart
cd $RUN_DIR

cat>$RUN_DIR/$SYS_NAME.lmp.inp<<EOF
###########################################################################
#
# Input file for uniaxial compressive loading
#
###########################################################################

# Restarting from a previus run $ID_RUN_RESTART
read_restart zca-ss-npt_PREV.restart

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

# set time integration
#
fix 1 all npt temp 300 300 1 y 0 0 1 z 0 0 1 drag 1

# set thermo output
thermo 1000
thermo_style custom step lx ly lz press pxx pyy pzz pe temp

# run one single step
run 1
unfix 1

# store final cell length for strain calculations
variable tmp equal "lx"
variable L0 equal \${tmp}
print "Initial Length, L0: \${L0}"

###########################################################################
#
# Deformation
#
###########################################################################

reset_timestep 0

fix 1 all npt temp 300 300 1 y 0 0 1 z 0 0 1 drag 1

# this is the strain rate in s^-1
variable srate equal $SRATE
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
fix def1 all print 100 "\${p1} \${p2} \${p3} \${p4}" file comp_100.def1.txt screen no

# display thermo
thermo $N_THERMO
thermo_style custom step v_strain temp v_p2 v_p3 v_p4 ke pe press vol lx ly lz

dump dss all atom $N_THERMO zca-ss.dump

# writing periodically a restart file
restart 10000000 zca-ss-npt.restart

run $N_RUN
EOF

$LMP_CMD < $RUN_DIR/$SYS_NAME.lmp.inp > $RUN_DIR/$SYS_NAME.lmp.out

rm $RUN_DIR/zca-ss-npt_PREV.restart

echo "Final do job:" `date`
echo "-----------------------------------------" 

