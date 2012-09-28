#!/bin/bash

# this script launches the specfem simulation
# Percy Galvez, ETH-zurich.Brutus.. 2010.

# use the normal queue unless otherwise directed

queue=""
if [ $# -eq 1 ]; then
	echo "Setting the queue to $1"
	queue="-q $1"
fi

d=`date`
echo "Starting compilation $d"
make clean
make

d=`date`

echo "Finished compilation $d"

# total number of nodes is the product of the values read
numnodes=`grep NPROC DATA/Par_file | cut -d = -f 2`

echo "Submitting job"
echo "bsub $queue -R ib  -n $numnodes -W 360 < go_specfem3D_lsf.bash"

# Intel version of open MPI needs infiniBand libraries (librdmacm.so),even
#therefore  nodes with this characteristic are required, in that sense option  
# -R ib will take infiniband nodes.
 
bsub -R ib -n $numnodes -W 360 < go_specfem3D_lsf.bash
