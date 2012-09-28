# Makefile

#############################################################
## Compiler : 
#F90 = ifort      # use -warn
F90 = gfortran    # use -Wall

## modify to match your library paths
SCOTCH_LIBS = -L/cluster/home/erdw/gpercy/scotch/scotch_5.1/lib/ -lscotch -lscotcherr
#(SCOTCH library address .Change it in case 
#your SCOTH is install in another address)

#############################################################

LIBS = fault_scotch.o\
               part_decompose_mesh_SCOTCH.o \
				decompose_mesh_SCOTCH.o \
				program_decompose_mesh_SCOTCH.o

# targets
all: xdecompose_mesh_SCOTCH

xdecompose_mesh_SCOTCH: $(LIBS)
	${F90} -o xdecompose_mesh_SCOTCH $(LIBS) $(SCOTCH_LIBS)

fault_scotch.o : fault_scotch.f90
	${F90} -c fault_scotch.f90

part_decompose_mesh_SCOTCH.o: part_decompose_mesh_SCOTCH.f90
	${F90} -c part_decompose_mesh_SCOTCH.f90

decompose_mesh_SCOTCH.o: decompose_mesh_SCOTCH.f90 part_decompose_mesh_SCOTCH.f90
	${F90} -c decompose_mesh_SCOTCH.f90

program_decompose_mesh_SCOTCH.o: program_decompose_mesh_SCOTCH.f90
	${F90} -c program_decompose_mesh_SCOTCH.f90

clean:
	rm -f *.o *.mod a.out xdecompose_mesh_SCOTCH
