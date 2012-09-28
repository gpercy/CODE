#!/usr/bin/env python
#    P. Galvez (ETH-Zurich, 25.10.2012): Tohoku mesh breaking the surface in a more realistic way.
import cubit
import cubit2specfem3d 
import math
import os
import sys
from save_fault_nodes_elements import *
from absorbing_boundary import *
from functions import m2km

cubit.cmd('reset')

### Reading mesh_japan.cub , meshed already in cubit ####.
### Remark , no need to open the fault , the fault already 
### has side up and down done with merge off ###
### WARNING : This mesh has been made in CUBIT 13.0 and 
### this mesh can not be read with previous CUBIT versions.

open_mesh = 'open "/Users/pgalvez/CODE/toho_real_mesh3.cub"' # The fault is closed. 
cubit.cmd(open_mesh)
cubit.cmd('compress ids hex face edge node')  
                                # THIS IS RELLY IMPORTANT TO AVOID GAPS IN NODES indexing  AND HAVING
                                # Indices bigger that the total grid points and 
                                # avoid fault segmentation fault errors
########### Slab elements and nodes ###############
########  SLAB  ################################################
os.system('mkdir -p MESH') 
Au = [358,366]  # A_up
Ad = [324,326]  # A_down

############ Fault opening ######################################
cubit.cmd('set node constrain off')                     # Here the mesh is in km. 
cubit.cmd('node in surf 358 366 move X 0 Y 0 Z 0.005') # delta= 10m = 10e-3 km.
cubit.cmd('node in surf 324 326 move X 0 Y 0 Z -0.005')  # the fault will be a delta in km in the z direction. 

##################################################################
# Saving fault elements.
faultA = fault_input(1,Au,Ad)
quads_Aup,quads_Adp = save_cracks(faultA.name,faultA.surface_u,faultA.surface_d)
#Unpacking list.
quads_Au=unpack_list(quads_Aup)
quads_Ad=unpack_list(quads_Adp)

print 'len(Au):',len(quads_Au)
print 'len(Ad):',len(quads_Ad)

if not (len(quads_Au)==len(quads_Ad)):
    print 'Number of elements for each fauld side up and down do not concide'
    sys.exit('goodbye')
   
save_elements_nodes(faultA.name,quads_Au,quads_Ad)

##  FOR THE BULK (Seismic wave propagation part for SESAME)
####### defining absorbing-boundary surface 
xmin   = [131,129,112,121,174,180]
xmax   = [406,407,408]
ymin   = [377,378,371,417,167,168,166,165]
ymax   = [137,136,125,124,376,389,391,412] 
bottom = [409,404,393,147,126,114,333,336,339,422,179,183,380,387,416]
#### defining free surface with block face_topo.
topo    = [175,190,373,414,385,420,343,357,367,116,411,402,397,143,134] # Free surface.
abs_surf = abs_surface_topo(xmin,xmax,ymin,ymax,bottom,topo)
##..... which extracts the bounding faces and defines them into blocks 
entities=['face']
## WARNING : The mesh is rotated respect to Z, 
##           the boundary surfaces are not parallel to cartesian coordinates.
define_bc_topo(entities,abs_surf) # Define absorbing boundaries surfaces and free surface. 
 
#### Define material properties for the volumes (18 volumes) ################ 
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################') 
 
# Material properties tohoku (Homogeneous). 
 
cubit.cmd('block 1 name "elastic 1" ')        # material region  
cubit.cmd('block 1 attribute count 5') 
cubit.cmd('block 1 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 1 attribute index 2 5800')   # vp 
cubit.cmd('block 1 attribute index 3 3420')    # vs 
cubit.cmd('block 1 attribute index 4 2700')   # rho 
cubit.cmd('block 1 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 2 name "elastic 2" ')        # material region  
cubit.cmd('block 2 attribute count 5') 
cubit.cmd('block 2 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 2 attribute index 2 5800')   # vp 
cubit.cmd('block 2 attribute index 3 3420')    # vs 
cubit.cmd('block 2 attribute index 4 2700')   # rho 
cubit.cmd('block 2 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 3 name "elastic 3" ')        # material region  
cubit.cmd('block 3 attribute count 5') 
cubit.cmd('block 3 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 3 attribute index 2 5800')   # vp 
cubit.cmd('block 3 attribute index 3 3420')    # vs 
cubit.cmd('block 3 attribute index 4 2700')   # rho 
cubit.cmd('block 3 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 4 name "elastic 4" ')        # material region  
cubit.cmd('block 4 attribute count 5') 
cubit.cmd('block 4 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 4 attribute index 2 5800')   # vp 
cubit.cmd('block 4 attribute index 3 3420')    # vs 
cubit.cmd('block 4 attribute index 4 2700')   # rho 
cubit.cmd('block 4 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 5 name "elastic 5" ')        # material region  
cubit.cmd('block 5 attribute count 5') 
cubit.cmd('block 5 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 5 attribute index 2 5800')   # vp 
cubit.cmd('block 5 attribute index 3 3420')    # vs 
cubit.cmd('block 5 attribute index 4 2700')   # rho 
cubit.cmd('block 5 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 6 name "elastic 6" ')        # material region  
cubit.cmd('block 6 attribute count 5') 
cubit.cmd('block 6 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 6 attribute index 2 5800')   # vp 
cubit.cmd('block 6 attribute index 3 3420')    # vs 
cubit.cmd('block 6 attribute index 4 2700')   # rho 
cubit.cmd('block 6 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 7 name "elastic 7" ')        # material region  
cubit.cmd('block 7 attribute count 5') 
cubit.cmd('block 7 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 7 attribute index 2 5800')   # vp 
cubit.cmd('block 7 attribute index 3 3420')    # vs 
cubit.cmd('block 7 attribute index 4 2700')   # rho 
cubit.cmd('block 7 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 8 name "elastic 8" ')        # material region  
cubit.cmd('block 8 attribute count 5') 
cubit.cmd('block 8 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 8 attribute index 2 5800')   # vp 
cubit.cmd('block 8 attribute index 3 3420')    # vs 
cubit.cmd('block 8 attribute index 4 2700')   # rho 
cubit.cmd('block 8 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 9 name "elastic 9" ')        # material region  
cubit.cmd('block 9 attribute count 5') 
cubit.cmd('block 9 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 9 attribute index 2 5800')   # vp 
cubit.cmd('block 9 attribute index 3 3420')    # vs 
cubit.cmd('block 9 attribute index 4 2700')   # rho 
cubit.cmd('block 9 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 10 name "elastic 10" ')        # material region  
cubit.cmd('block 10 attribute count 5') 
cubit.cmd('block 10 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 10 attribute index 2 5800')   # vp 
cubit.cmd('block 10 attribute index 3 3420')    # vs 
cubit.cmd('block 10 attribute index 4 2700')   # rho 
cubit.cmd('block 10 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 11 name "elastic 11" ')        # material region  
cubit.cmd('block 11 attribute count 5') 
cubit.cmd('block 11 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 11 attribute index 2 5800')   # vp 
cubit.cmd('block 11 attribute index 3 3420')    # vs 
cubit.cmd('block 11 attribute index 4 2700')   # rho 
cubit.cmd('block 11 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 12 name "elastic 12" ')        # material region  
cubit.cmd('block 12 attribute count 5') 
cubit.cmd('block 12 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 12 attribute index 2 5800')   # vp 
cubit.cmd('block 12 attribute index 3 3420')    # vs 
cubit.cmd('block 12 attribute index 4 2700')   # rho 
cubit.cmd('block 12 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 13 name "elastic 13" ')        # material region  
cubit.cmd('block 13 attribute count 5') 
cubit.cmd('block 13 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 13 attribute index 2 5800')   # vp 
cubit.cmd('block 13 attribute index 3 3420')    # vs 
cubit.cmd('block 13 attribute index 4 2700')   # rho 
cubit.cmd('block 13 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 14 name "elastic 14" ')        # material region  
cubit.cmd('block 14 attribute count 5') 
cubit.cmd('block 14 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 14 attribute index 2 5800')   # vp 
cubit.cmd('block 14 attribute index 3 3420')    # vs 
cubit.cmd('block 14 attribute index 4 2700')   # rho 
cubit.cmd('block 14 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 15 name "elastic 15" ')        # material region  
cubit.cmd('block 15 attribute count 5') 
cubit.cmd('block 15 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 15 attribute index 2 5800')   # vp 
cubit.cmd('block 15 attribute index 3 3420')    # vs 
cubit.cmd('block 15 attribute index 4 2700')   # rho 
cubit.cmd('block 15 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 16 name "elastic 16" ')        # material region  
cubit.cmd('block 16 attribute count 5') 
cubit.cmd('block 16 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 16 attribute index 2 5800')   # vp 
cubit.cmd('block 16 attribute index 3 3420')    # vs 
cubit.cmd('block 16 attribute index 4 2700')   # rho 
cubit.cmd('block 16 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 17 name "elastic 17" ')        # material region  
cubit.cmd('block 17 attribute count 5') 
cubit.cmd('block 17 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 17 attribute index 2 5800')   # vp 
cubit.cmd('block 17 attribute index 3 3420')    # vs 
cubit.cmd('block 17 attribute index 4 2700')   # rho 
cubit.cmd('block 17 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 18 name "elastic 18" ')        # material region  
cubit.cmd('block 18 attribute count 5') 
cubit.cmd('block 18 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 18 attribute index 2 5800')   # vp 
cubit.cmd('block 18 attribute index 3 3420')    # vs 
cubit.cmd('block 18 attribute index 4 2700')   # rho 
cubit.cmd('block 18 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 19 name "elastic 19" ')        # material region  
cubit.cmd('block 19 attribute count 5') 
cubit.cmd('block 19 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 19 attribute index 2 5800')   # vp 
cubit.cmd('block 19 attribute index 3 3420')    # vs 
cubit.cmd('block 19 attribute index 4 2700')   # rho 
cubit.cmd('block 19 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 20 name "elastic 20" ')        # material region  
cubit.cmd('block 20 attribute count 5') 
cubit.cmd('block 20 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 20 attribute index 2 5800')   # vp 
cubit.cmd('block 20 attribute index 3 3420')    # vs 
cubit.cmd('block 20 attribute index 4 2700')   # rho 
cubit.cmd('block 20 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 21 name "elastic 21" ')        # material region  
cubit.cmd('block 21 attribute count 5') 
cubit.cmd('block 21 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 21 attribute index 2 5800')   # vp 
cubit.cmd('block 21 attribute index 3 3420')    # vs 
cubit.cmd('block 21 attribute index 4 2700')   # rho 
cubit.cmd('block 21 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 22 name "elastic 22" ')        # material region  
cubit.cmd('block 22 attribute count 5') 
cubit.cmd('block 22 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 22 attribute index 2 5800')   # vp 
cubit.cmd('block 22 attribute index 3 3420')    # vs 
cubit.cmd('block 22 attribute index 4 2700')   # rho 
cubit.cmd('block 22 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 23 name "elastic 23" ')        # material region  
cubit.cmd('block 23 attribute count 5') 
cubit.cmd('block 23 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 23 attribute index 2 5800')   # vp 
cubit.cmd('block 23 attribute index 3 3420')    # vs 
cubit.cmd('block 23 attribute index 4 2700')   # rho 
cubit.cmd('block 23 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )


#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH') 

## transforming coordinates from km to m .
m2km()

# Fault opening has been done in CUBIT.
# INPUTS :
#file_nodes_coord = 'MESH/nodes_coords_file'
#fault_file =       'MESH/fault_file_1.dat'
# OUTPUTS :
#name_out = 'MESH/nodes_coords_file_open_fault'
#fsideu   = 'MESH/fault_sideu.dat'
#fsided   = 'MESH/fault_sided.dat'
#delta = 2.0 # Make sure that this delta is in coorcondance with 
            # FAULT_GAP_TOLERANCE which is normally set up to 1.0d0 but can be changed 
            # in decompose_mesh_SCOTH/fault_scotch.f90, and higher
            # than constants.h/SMALLVAL_TOL used by the routine get_global.f90.
            # SMALLVALTOL = 1.d-10 * dabs(UTM_X_MAX - UTM_X_MIN)

#nodes_coords_fault_open(file_nodes_coord,fault_file,name_out,fsideu,fsided,delta)

# all files needed by SCOTCH are now in directory MESH 

