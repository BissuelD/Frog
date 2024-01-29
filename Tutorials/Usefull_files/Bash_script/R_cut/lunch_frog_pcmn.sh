#!/bin/bash
#
### variables SGE
#
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N gamma_frog
### file d'attente (a changer)

#$ -q  h48*,h6*
#$ -pe openmp8 8

### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V
### mails en debut et fin d'execution
# #$ -m be
# aller dans le repertoire de travail/soumission
# important, sinon, le programme est lancé depuis ~/

#NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')
#echo "Total CPU count = $NP"

cd ${SGE_O_WORKDIR}
### configurer l'environnement (a changer)
module purge
module load Python/3.6.1
source /home/glebreto/Software/Frog/vfrog/bin/activate

# export LD_LIBRARY_PATH=/applis/PSMN/debian9/software/Core/Python/3.6.1/bin/python:$LD_LIBRARY_PATH:/applis/PSMN/debian9/software/Core/Python/3.6.1/bin/python 

# MPIRUN="/usr/bin/mpirun"
MPIRUN="/usr/bin/mpirun.openmpi"

echo ${MPIRUN}
echo ${TMPDIR}


HOSTFILE=${TMPDIR}/machines

echo ${MPIRUN}

echo `which python`
# export OMP_NUM_THREADS=1

# $MPIRUN -v -hostfile ${HOSTFILE} -np 1 Frog parameters_frog.py > frog_output_${JOB_ID} # for mpirun
$MPIRUN -np 1 Frog parameters_frog.py > frog_output_${JOB_ID}

# $MPIRUN -np 1 Frog parameters_frog.py > frog_output_${JOB_ID}
# "${MPIRUN}" -v -x LD_LIBRARY_PATH -hostfile "${HOSTFILE}" -np 1 Frog parameters_frog.py > frog_output_${JOB_ID}


