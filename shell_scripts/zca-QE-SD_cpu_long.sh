#!/bin/bash
#SBATCH --nodes=8                   # Numero de Nos
#SBATCH --ntasks-per-node=24        # Numero de tarefas por No
#SBATCH --ntasks=192                # Numero total de tarefas MPI
#SBATCH -p cpu_long                # Fila (partition) a ser utilizada
#SBATCH -J QE-47476                # Nome job
#SBATCH --exclusive                 # Utilizacao exclusiva dos nos durante a execucao do job
#SBATCH -o stdout.%j                # File to which STDOUT will be written
#SBATCH -e stderr.%j                # File to which STDERR will be written
#SBATCH --mail-type=END             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent
#SBATCH --time=10:00:00             # Altera o tempo limite HH:MM:SS

# Exibe os nos alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

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
module load mathlibs/fftw/3.3.8_intel
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

#
# Command lines
#
export CMD_QE="srun -n 192 $SCRATCH/opt/qe-6.2_ARY_FFTW3_EXT_INTEL/bin/pw.x -npool 8"

### Variaveis de controle da simulacao 1 ###
export NAT=100
export BD_DIR="MG-NMR/ML/big-data-full"
export SYS_NAME="${1}"

export ID_RUN="${2}"
export SUB_STEP="${3}"

export BASE_DIR="$BD_DIR/$SYS_NAME/c/md/lammps/$NAT/$ID_RUN/2000"
export WORK_DIR="$SCRATCH/work/$BASE_DIR/$SUB_STEP"
export RUN_DIR="$SCRATCH/$BASE_DIR/$SUB_STEP"

cd $RUN_DIR

echo $QE_RUN > $RUN_DIR/QE_run

rm -rf $WORK_DIR
mkdir -p $WORK_DIR

$CMD_QE < $RUN_DIR/$SYS_NAME.scf.in > $RUN_DIR/$SYS_NAME.scf.out

### Variaveis de controle da simulacao 2 ###
export NAT=100
export BD_DIR="MG-NMR/ML/big-data-full"
export SYS_NAME="${4}"

export ID_RUN="${5}"
export SUB_STEP="${6}"

export BASE_DIR="$BD_DIR/$SYS_NAME/c/md/lammps/$NAT/$ID_RUN/2000"
export WORK_DIR="$SCRATCH/work/$BASE_DIR/$SUB_STEP"
export RUN_DIR="$SCRATCH/$BASE_DIR/$SUB_STEP"

cd $RUN_DIR

echo $QE_RUN > $RUN_DIR/QE_run

rm -rf $WORK_DIR
mkdir -p $WORK_DIR

$CMD_QE < $RUN_DIR/$SYS_NAME.scf.in > $RUN_DIR/$SYS_NAME.scf.out

