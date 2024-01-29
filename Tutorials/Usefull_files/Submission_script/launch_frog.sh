#!/bin/bash
#$ -S /bin/bash
#$ -N yourname
#$ -q  h48*
#$ -pe openmp8 8

#NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')
#echo "Total CPU count = $NP"

cd ${SGE_O_WORKDIR}
### configurer l'environnement (a changer)
module purge
module load Python/3.6.1
source /home/glebreto/Software/Frog/vfrog/bin/activate

MPIRUN="/usr/bin/mpirun.openmpi"

echo ${MPIRUN}
echo ${TMPDIR}

HOSTFILE=${TMPDIR}/machines

echo ${MPIRUN}

echo `which python`

# 8 workers are available here: the parralellization is set in Frog input file and not according to the below command line: 
$MPIRUN -np 1 Frog parameters_frog.py > frog_output_${JOB_ID}


