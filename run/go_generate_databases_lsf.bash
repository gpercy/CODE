#!/bin/bash

## job name and output file
#BSUB -J go_generate_databases
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -e OUTPUT_FILES/%J.e

###########################################################

if [ -z $USER ]; then
	echo "could not run go_mesher_...bash as no USER env is set"
	exit 2
fi

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

cp DATA/Par_file OUTPUT_FILES/

# obtain lsf job information
cat $LSB_DJOB_HOSTFILE > OUTPUT_FILES/compute_nodes
echo "$LSB_JOBID" > OUTPUT_FILES/jobid

echo starting MPI mesher on $numnodes processors
echo " "

sleep 2 
mpirun $PWD/xgenerate_databases

echo "done "
