#!/bin/bash
make clean
make 
# ./xdecompose_mesh_SCOTCH 4 MESH/OUTPUT_FILES ../DATABASES_MPI/
echo "./xdecompose_mesh_SCOTCH $npart $input_dir $ouput_dir "
./xdecompose_mesh_SCOTCH 100 ~/FAULT_SOURCE/EXAMPLES/splay_faults/MESH /cluster/work/erdw/gpercy
