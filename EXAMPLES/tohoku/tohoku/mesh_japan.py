#!/usr/bin/env python

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

open_mesh = 'open "/Users/pgalvez/TOHOKU/MESH/mesh_japan_ok.cub"'
cubit.cmd(open_mesh)

cubit.cmd('compress ids hex face edge node')  
                                # THIS IS RELLY IMPORTANT TO AVOID GAPS IN NODES indexing  AND HAVING
                                # Indices bigger that the total grid points and 
                                # avoid fault segmentation fault errors

########### Slab elements and nodes ###############
########  SLAB  ################################################
os.system('mkdir -p MESH') 

Au = [42]  # A_up
Ad = [36]  # A_down

####### Fault opening #############################################
cubit.cmd('set node constraint off')
cubit.cmd('node in surf 42 move X 0 Y 0 Z 0.001')   # delta = 2e-3 km. in Z direction.
cubit.cmd('node in surf 36 move X 0 Y 0 Z -0.001')  # In general the shift should be normal to the fault.

####################################################################
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
####### defining absorbing-boundary surface...
xmin = [133,114,124,97,129,106]
xmax = [175,151,170,144,180,161]
ymin = [79,105,128,53,154,176]
ymax = [179,87,135,117,63,160] 
abs_surf = abs_surface(xmin,xmax,ymin,ymax)
## Fixing the bottom to find bottom absorbing boundaries.
zmin = -150
zmax = 0 
##..... which extracts the bounding faces and defines them into blocks 
entities=['face']
## WARNING : The mesh is rotated respect to Z, 
##           the boundary surfaces are not parallel to cartesian coordinates.
define_bc(entities,zmin,zmax,abs_surf) 
 
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

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH') 

m2km()
### this mesh is in km . It is more convinient to work in meters 
### changing nodes_coords_file to meters . 


 
 
# all files needed by SCOTCH are now in directory MESH 
