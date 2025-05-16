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
#SBATCH --time=10:00:00             # Altera o tempo limite HH:MM:SS

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
export LMP_CMD="srun -n 504 /scratch/rmnvm2/ary.junior2/opt/lammps-16Feb16_ARY/src/lmp_mpi_fftw3_SD_update"

### Variaveis de controle da simulacao ###
export SYS_NAME="${1}"
export BASE_DIR="/scratch/rmnvm2/ary.junior2/MG-NMR/SS-ML"

### Rodando a simulacao ###
export GDIR="${2}"
#export GDIR="`cat $BASE_DIR/GDIR`"
export RUN_DIR="$BASE_DIR/$SYS_NAME/c/md/lammps/ultimate/2/$GDIR"
mkdir -p $RUN_DIR

cp $BASE_DIR/ZrCuAl-2011.lammps.eam $RUN_DIR
cp $BASE_DIR/$SYS_NAME/$SYS_NAME-10000-1-1-1.data $RUN_DIR
cd $RUN_DIR

$LMP_CMD < $RUN_DIR/$SYS_NAME.lmp.inp > $RUN_DIR/$SYS_NAME.lmp.out

echo "Final do job:" `date`
echo "-----------------------------------------" 

