#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d 
import math
import os
import sys
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

fault_u = 2  # fault surface up   : surface 2
fault_d = 3  # fault surface down : surface 3
txt =''

fault_file=open('MESH/fault_file_1.dat','w')
list_hex=cubit.parse_cubit_list('hex','all')

quads_fault_u = cubit.get_surface_quads(fault_u)
quads_fault_d = cubit.get_surface_quads(fault_d)

# TO DO : stop python properly in case fault nodes at both sides
#         do not match.
if len(quads_fault_u) != len(quads_fault_d):
   stop
# Writting number of elements at both sides to make 
# double sure everything is going fine .
txt='%10i %10i\n' % (len(quads_fault_u),len(quads_fault_d))
fault_file.write(txt)

dic_quads_fault_u = dict(zip(quads_fault_u,quads_fault_u)) 
dic_quads_fault_d = dict(zip(quads_fault_d,quads_fault_d)) 


# FAULT SIDE DOWN
for h in list_hex: 
    faces = cubit.get_sub_elements('hex',h,2)  
    for f in faces:
        if dic_quads_fault_d.has_key(f): 
           nodes=cubit.get_connectivity('Face',f)
#           print 'h,fault nodes side down :',h,nodes[0],nodes[1],nodes[2],nodes[3]
           txt='%10i %10i %10i %10i %10i\n' % (h,nodes[0],\
                                             nodes[1],nodes[2],nodes[3])
           fault_file.write(txt)

# FAULT SIDE UP
for h in list_hex: 
    faces = cubit.get_sub_elements('hex',h,2)  
    for f in faces:
        if dic_quads_fault_u.has_key(f): 
           nodes=cubit.get_connectivity('Face',f)
#           print 'h,fault nodes side up :',h,nodes[0],nodes[1],nodes[2],nodes[3]
           txt='%10i %10i %10i %10i %10i\n' % (h,nodes[0],\
                                             nodes[1],nodes[2],nodes[3])
           fault_file.write(txt)

# CLOSING FAULT FILE 
fault_file.close()

#  FOR THE BULK (Seismic wave propagation part for SESAME)

###### This is boundary_definition.py of GEOCUBIT 
#..... which extracts the bounding faces and defines them into blocks 
boundary_definition.entities=['face'] 
boundary_definition.define_bc(boundary_definition.entities,parallel=True) 
 
#### Define material properties for the 2 volumes ################ 
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

