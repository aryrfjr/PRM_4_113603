#!/bin/bash
#SBATCH --nodes=4                   # Numero de Nos
#SBATCH --ntasks-per-node=24        # Numero de tarefas por No
#SBATCH --ntasks=96                # Numero total de tarefas MPI
#SBATCH -p cpu_dev                # Fila (partition) a ser utilizada
#SBATCH -J zca-npt                # Nome job
#SBATCH --exclusive                 # Utilizacao exclusiva dos nos durante a execucao do job
#SBATCH -o stdout.%j                # File to which STDOUT will be written
#SBATCH -e stderr.%j                # File to which STDERR will be written
#SBATCH --mail-type=END             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent
#SBATCH --time=00:20:00             # Altera o tempo limite HH:MM:SS

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
export LMP_CMD="srun -n 96 /scratch/rmnvm/ary.junior/opt/lammps-16Feb16_ARY/src/lmp_mpi_fftw3_SD_update"

### Variaveis de controle da simulacao ###
export ID_RUN=${1}
export SYS_NAME="${2}"
export BASE_DIR="/scratch/rmnvm/ary.junior/MG-NMR/SS-ML"
export ID_RUN_RESTART=${3}

### Rodando a simulacao ###
export RUN_DIR="$BASE_DIR/$SYS_NAME/c/md/lammps/$ID_RUN"
mkdir -p $RUN_DIR

cp $BASE_DIR/ZrCuAl-2011.lammps.eam $RUN_DIR
cp $BASE_DIR/$SYS_NAME/c/md/lammps/$ID_RUN_RESTART/zca-equ-nvt-300.restart $RUN_DIR
cd $RUN_DIR

cat>$RUN_DIR/$SYS_NAME.lmp.inp<<EOF
###########################################################################
#
# Equilibration of the loaded atoms in NPT
#
###########################################################################

# Restarting from the structure in equilibrium at 300 K from run $ID_RUN_RESTART
read_restart zca-equ-nvt-300.restart

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

#
# And then the NPT ensemble to equilibrate pressure.

# setting time integration
#
fix 1 all npt temp 300 300 1 iso 0 0 1 drag 1

# set thermo output
thermo 1000
thermo_style custom step lx ly lz press pxx pyy pzz pe temp

# run for 100 picoseconds
run 100000

write_restart zca-equ-npt-300.restart
EOF

$LMP_CMD < $RUN_DIR/$SYS_NAME.lmp.inp > $RUN_DIR/$SYS_NAME.lmp.out

rm $RUN_DIR/zca-equ-nvt-300.restart

echo "Final do job:" `date`
echo "-----------------------------------------" 

