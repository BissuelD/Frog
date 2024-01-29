#!/bin/bash
#
### variables SGE
#
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N wat_short_run_approx_15min
### file d'attente (a changer)
#$ -q CLG52*,CLG6*,r815lin128ib,E5-2697Av4deb256*,SLG6142* 
# $ -q x5*
# $ -q h6-E5-2667v4deb128,h48-E5-2670deb128 
### parallel environnement & nslots (a changer)
#$ -pe mpi32_debian 32
# $ -pe mpi16_debian 16  
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
module load Dalton/2018.2-full-debian9
export OMP_NUM_THREADS=1

MPIRUN="/usr/bin/mpirun"
HOSTFILE=${TMPDIR}/machines



# Unix commande to perform the QM runs 

SCRATCH_DIR=/tmp/glebreto


