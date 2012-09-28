#!/usr/bin/env python

import cubit
import cubit2specfem3d 
import math
import os
import sys
from save_fault_nodes_elements import *
from absorbing_boundary import *

cubit.cmd('reset')

km = 1000
z_surf=0*km
####  initializing coordinates x,y,z
x=[]     # fault
y=[]
z=[]

xbulk=[] # bulk
ybulk=[]
zbulk=[]

xbulk.append(-21*km)   #x1
xbulk.append(21*km)    #x2
xbulk.append(21*km)    #x3
xbulk.append(-21*km)   #x4

ybulk.append(-21*km)  #y1
ybulk.append(-21*km)  #y2
ybulk.append(21*km)   #y3
ybulk.append(21*km)   #y4

zbulk=[z_surf]*4

x.append(-9*km) #x5
x.append(0*km)   #x6
x.append(9*km)  #x7
x.append(0*km)   #x8

y.append(0.0)       #y5
y.append(0.1)    #y6
y.append(0.0)       #y7
y.append(-0.1)   #y8


z=[z_surf]*4

####################  bulk ###########################################
for i in range(len(xbulk)): 
   vert="create vertex x "+str(xbulk[i])+" y "+str(ybulk[i])+" z "+str(zbulk[i]) 
   cubit.cmd(vert) 

################  Loading fault points profile#############################
for i in range(len(x)):
  vert="create vertex x "+str(x[i])+" y "+str(y[i])+" z "+str(z[i])
  cubit.cmd(vert)


################ creating fault domains #################################
bulk1="create curve vertex 1 2"
bulk2="create curve vertex 2 3"
bulk3="create curve vertex 3 4"
bulk4="create curve vertex 4 1"
 
fault_up="create curve spline vertex 5 6 7"
fault_down="create curve spline vertex 5 8 7" 

cubit.cmd(bulk1)
cubit.cmd(bulk2)
cubit.cmd(bulk3)
cubit.cmd(bulk4)

cubit.cmd(fault_up) 
cubit.cmd(fault_down) 

surface="create surface curve 1 2 3 4 5 6"
cubit.cmd(surface)

cubit.cmd("sweep surface 1  vector 0  0 -1 distance "+str(21*km)) 
cubit.cmd("sweep curve 5 vector 0 0 -1 distance "+str(21*km)) 
cubit.cmd("sweep curve 6 vector 0 0 -1 distance "+str(21*km)) 

#####################################################
elementsize = 1000
cubit.cmd("imprint all")
cubit.cmd("merge all")
cubit.cmd("surface 1 size "+str(elementsize))
cubit.cmd("volume 1 size "+str(elementsize))
cubit.cmd("surface 1 scheme pave")
cubit.cmd("mesh surface 1")
cubit.cmd("mesh volume 1")
#cubit.cmd("unmerge surface 2 3")

########### SIDESETS (NOT USED ) ###############
#fault_A_elements_up="sideset 200 surface 2"
#fault_A_elements_down="sideset 201 surface 3"
#cubit.cmd(fault_A_elements_up)
#cubit.cmd(fault_A_elements_down)

########### Fault elements and nodes ###############

os.system('mkdir -p MESH') 

Au = [2]  # A_up
Ad = [3]  # A_down

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
xmin = [7]
xmax = [5]
ymin = [4]
ymax = [6] 
abs_surf = abs_surface(xmin,xmax,ymin,ymax)
## Fixing the bottom to find bottom absorbing boundaries.
zmin = -21000
zmax = 0 
##..... which extracts the bounding faces and defines them into blocks 
entities=['face']
## WARNING : The mesh is rotated respect to Z, 
##           the boundary surfaces are not parallel to cartesian coordinates.
define_bc(entities,zmin,zmax,abs_surf) 
 
#### Define material properties for the volumes (18 volumes) ################ 
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################') 

# Material properties in concordance with tpv5 benchmark. 
 
cubit.cmd('block 1 name "elastic 1" ')        # material region  
cubit.cmd('block 1 attribute count 5') 
cubit.cmd('block 1 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 1 attribute index 2 6000')   # vp 
cubit.cmd('block 1 attribute index 3 3464')    # vs 
cubit.cmd('block 1 attribute index 4 2670')   # rho 
cubit.cmd('block 1 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH')  
 
# all files needed by SCOTCH are now in directory MESH 

