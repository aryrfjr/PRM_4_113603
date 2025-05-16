#!/bin/bash
#SBATCH --nodes=4                   # Numero de Nos
#SBATCH --ntasks-per-node=24        # Numero de tarefas por No
#SBATCH --ntasks=96                # Numero total de tarefas MPI
#SBATCH -p cpu_long                # Fila (partition) a ser utilizada
#SBATCH -J zca-46468                # Nome job
#SBATCH --exclusive                 # Utilizacao exclusiva dos nos durante a execucao do job
#SBATCH -o stdout.%j                # File to which STDOUT will be written
#SBATCH -e stderr.%j                # File to which STDERR will be written
#SBATCH --mail-type=END             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ary.ferreira@df.ufscar.br  # Email to which notifications will be sent
#SBATCH --time=300:00:00             # Altera o tempo limite HH:MM:SS

#Exibe os nós alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

echo "-----------------------------------------"
echo "Inicio do job:" `date`

#
# Modulos
#
source /scratch/app/modulos/intel-psxe-2017.sh
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/app/fftw/3.3.7_intel/lib

### Definir variaveis para processamento com MPI da SGI ###

export MPI_MEMMAP_OFF=1
export MPI_GROUP_MAX=1000
export MPI_COMM_MAX=1000
export MPI_TYPE_MAX=100000

### A opcao -np indica o numero total de processos.  ###
export LMP_CMD="srun -n 96 /scratch/rmnvm/ary.junior/opt/lammps-16Feb16_ARY/src/lmp_mpi_fftw3"

### Variaveis de controle da simulacao ###
export SYMBS=(Zr Cu Al)
export NATS=(46 46 8)
export NAT=100
export BD_DIR="/scratch/rmnvm/ary.junior/MG-NMR/ML/big-data-full"
export SYS_NAME="${SYMBS[0]}${NATS[0]}${SYMBS[1]}${NATS[1]}${SYMBS[2]}${NATS[2]}"
export LID_RUN=`cat $BD_DIR/last_id_run-$SYS_NAME`
export NRUNS=`cat $BD_DIR/n_runs-$SYS_NAME`
export BASE_DIR="$BD_DIR/$SYS_NAME/c/md/lammps/$NAT"
export CELL=6.98

### Rodando a simulacao ###
for ID_RUN in `seq $((LID_RUN+1)) $((LID_RUN+NRUNS))`
do

# Numero randomico utilizado pelo LAMMPS
export RNDATS=($(shuf -i 0-100000 -n 3))

export RUN_DIR="$BASE_DIR/$ID_RUN"
mkdir -p $RUN_DIR

cp $BD_DIR/ZrCuAl-2011.lammps.eam $RUN_DIR
cd $RUN_DIR

cat>$RUN_DIR/$SYS_NAME.lmp.inp<<EOF
# General setup
units metal
boundary p p p
atom_style atomic

# Initial random distribution of $NAT atoms
region box block -$CELL $CELL -$CELL $CELL -$CELL $CELL 
create_box 3 box
create_atoms 1 random ${NATS[0]} ${RNDATS[0]} NULL units box # ${SYMBS[0]}
create_atoms 2 random ${NATS[1]} ${RNDATS[1]} NULL units box # ${SYMBS[1]}
create_atoms 3 random ${NATS[2]} ${RNDATS[2]} NULL units box # ${SYMBS[2]}

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

# Summary of thermodynamic info
thermo 50
thermo_style custom step etotal temp press pe vol density lx ly lz enthalpy

# CG minimization at constant volume
# just to avoid overlapping of atoms
min_style cg # Conjugate gradient
minimize 0.0 1.0e-8 10000 10000

# Now trying to reach a liquid in thermal
# equilibrium at 2000 K.
# See more on 03/11/2016-(10).
fix th all npt temp 2000 2000 0.2 iso 0 0 2
# A time step of 2 fs.
timestep 0.002
thermo 5000
run 1000000
unfix th

#
# Now the quenching at 8.5e9 K/s from 2000 K to 300 K
#
fix thq all npt temp 2000 300 0.2 iso 0 0 2
thermo 50000
run 100000000
unfix thq

#
# Now the second thermalization at 300 K
#
fix th300 all npt temp 300 300 0.2 iso 0 0 2
dump zca-th300 all atom 500 zca-th300.dump
thermo 5000
run 1000000
write_restart zca-th300.restart
unfix th300
undump zca-th300
EOF

$LMP_CMD < $RUN_DIR/$SYS_NAME.lmp.inp > $RUN_DIR/$SYS_NAME.lmp.out

echo $ID_RUN > $BD_DIR/last_id_run-$SYS_NAME

done

echo "Final do job:" `date`
echo "-----------------------------------------" 

