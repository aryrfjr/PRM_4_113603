#!/bin/bash
#SBATCH --nodes=21                   # Numero de Nos
#SBATCH --ntasks-per-node=24        # Numero de tarefas por No
#SBATCH --ntasks=504                # Numero total de tarefas MPI
#SBATCH -p cpu                # Fila (partition) a ser utilizada
#SBATCH -J LOB-zca                # Nome job
#SBATCH --exclusive                 # Utilizacao exclusiva dos nos durante a execucao do job
#SBATCH -o stdout.%j                # File to which STDOUT will be written
#SBATCH -e stderr.%j                # File to which STDERR will be written
#SBATCH --mail-type=END             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent
#SBATCH --time=48:00:00             # Altera o tempo limite HH:MM:SS

# Exibe os nos alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

#
# MPI and Multithreads
#
export OMP_NUM_THREADS=24
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
export CMD_LB="srun --resv-ports --nodes 1 --ntasks=1 $SCRATCH/opt/lobster-3.0.0/lobster-3.0.0"

### Variaveis de controle da simulacao ###
export NAT=100
export BD_DIR="MG-NMR/ML/big-data-full"

export SYS_NAME="${1}"
export LOB_RUN="${2}"

while IFS=" " read -r ID_RUN SUB_STEP remainder
do

export BASE_DIR="$BD_DIR/$SYS_NAME/c/md/lammps/$NAT/$ID_RUN/2000"
export WORK_DIR="$SCRATCH/work/$BASE_DIR/$SUB_STEP"
export RUN_DIR="$SCRATCH/$BASE_DIR/$SUB_STEP"

cd $RUN_DIR

echo $LOB_RUN > $RUN_DIR/LOB_run

$CMD_LB > $RUN_DIR/$SYS_NAME.lb.out &

done < $SCRATCH/$BD_DIR/LOB_RUNs/lobr-$SYS_NAME-cpu-$LOB_RUN

# http://www.bosontreinamentos.com.br/linux/certificacao-lpic-1/comandos-wait-e-sleep-temporizacao-de-comandos-no-linux
wait

# TODO: loop again to remove
#rm -rf $WORK_DIR

