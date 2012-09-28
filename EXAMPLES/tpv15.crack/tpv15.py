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
R = math.pi/180
q = math.sqrt(3)
h = 0.1 

L = 64*km
####  initializing coordinates x,y,z
x=[]     # fault
y=[]
z=[]

xbulk=[] # bulk
ybulk=[]
zbulk=[]

xbulk.append(-L/2)   #x1
xbulk.append(L/2)    #x2
xbulk.append(L/2)    #x3
xbulk.append(-L/2)   #x4

ybulk.append(-L/2)  #y1
ybulk.append(-L/2)  #y2
ybulk.append(L/2)   #y3
ybulk.append(L/2)   #y4

zbulk=[z_surf]*4

### Main Fault #######################

x.append(-16*km) #x5
x.append(-50/math.tan(R*30))   #x6
x.append(L/2)  #x7
x.append(-50/math.tan(R*30))   #x8

y.append(0.0)    #y5
y.append(h)      #y6
y.append(0.0)    #y7
y.append(-h)     #y8

### Branch Fault ######################

x.append(x[3])                        #x9 = x8  (triple joint)
x.append((L/2/q)*math.cos(R*30))         #x10
x.append((L/q)*math.cos(R*30))        #x11
x.append((L/2/q)*math.cos(R*30))         #x12

y.append(y[3])                         #y9 = y8 (triple joint)
y.append(-50-(L/2/q)*math.sin(R*30)+h)   #y10
y.append(-50-(L/q)*math.sin(R*30))    #y11
y.append(-50-(L/2/q)*math.sin(R*30)-h)   #y12

z=[z_surf]*12

####################  bulk ###########################################
for i in range(len(xbulk)): 
   vert="create vertex x "+str(xbulk[i])+" y "+str(ybulk[i])+" z "+str(zbulk[i]) 
   cubit.cmd(vert) 

################  Loading fault points profile#############################
for i in range(len(x)):
  vert="create vertex x "+str(x[i])+" y "+str(y[i])+" z "+str(z[i])
  cubit.cmd(vert)


################ creating fault domains #################################
bulk1="create curve vertex 1 2"   # c1
bulk2="create curve vertex 2 11"  # c2
bulk3="create curve vertex 11 7"  # c3
bulk4="create curve vertex 7 3"   # c4
bulk5="create curve vertex 3 4"   # c5
bulk6="create curve vertex 4 1"   # c6

cubit.cmd(bulk1)
cubit.cmd(bulk2)
cubit.cmd(bulk3)
cubit.cmd(bulk4)
cubit.cmd(bulk5)
cubit.cmd(bulk6)

#### Main Fault ############################## 
fault_up_r="create curve spline vertex 5 6 "    #c7
fault_up_l="create curve spline vertex 6 7"     #c8

fault_down_r="create curve spline vertex 5 8"   #c9 
fault_down_l="create curve spline vertex 8 7"   #c10


cubit.cmd(fault_up_r) 
cubit.cmd(fault_up_l) 
cubit.cmd(fault_down_r)
cubit.cmd(fault_down_l)
 

#### Branch Fault ############################# 
fault_up="create curve spline vertex 9 10 11"    #c11 
fault_down="create curve spline vertex 9 12 11"  #c12

cubit.cmd(fault_up) 
cubit.cmd(fault_down) 
###############################################

surface1="create surface curve 2 1 6 5 4 7 8 9 12"

cubit.cmd(surface1)
surface2="create surface curve 3 10 11"
cubit.cmd(surface2)

cubit.cmd("sweep surface 1  vector 0  0 -1 distance "+str(21*km)) 
cubit.cmd("sweep surface 2  vector 0  0 -1 distance "+str(21*km)) 

#cubit.cmd("sweep curve 5 vector 0 0 -1 distance "+str(21*km)) 
#cubit.cmd("sweep curve 6 vector 0 0 -1 distance "+str(21*km)) 

#####################################################
elementsize = 300
cubit.cmd("imprint all")
cubit.cmd("merge all")
cubit.cmd("surface 1 2 size "+str(elementsize))
cubit.cmd("volume 1 2 size "+str(elementsize))
cubit.cmd("surface 1 2 scheme pave")
cubit.cmd("mesh surface 1 2")
cubit.cmd("mesh volume 1 2")

########### Fault elements and nodes ###############

## Main Fault ####################################################################
os.system('mkdir -p MESH') 

Au = [6,7]   # A_up
Ad = [5,14]  # A_down

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

####### Branch Fault##############################################################################
Bu = [4]   # B_up
Bd = [15]  # B_down

faultB = fault_input(2,Bu,Bd)
quads_Bup,quads_Bdp = save_cracks(faultB.name,faultB.surface_u,faultB.surface_d)
#Unpacking list.
quads_Bu=unpack_list(quads_Bup)
quads_Bd=unpack_list(quads_Bdp)

print 'len(Au):',len(quads_Bu)
print 'len(Ad):',len(quads_Bd)

if not (len(quads_Bu)==len(quads_Bd)):
    print 'Number of elements for each fauld side up and down do not concide'
    sys.exit('goodbye')
   
save_elements_nodes(faultB.name,quads_Bu,quads_Bd)

##  FOR THE BULK (Seismic wave propagation part for SESAME)

####### This is boundary_definition.py of GEOCUBIT 
##..... which extracts the bounding faces and defines them into blocks 
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


# Material properties in concordance with tpv5 benchmark. 
 
cubit.cmd('block 2 name "elastic 2" ')        # material region  
cubit.cmd('block 2 attribute count 5') 
cubit.cmd('block 2 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 2 attribute index 2 6000')   # vp 
cubit.cmd('block 2 attribute index 3 3464')    # vs 
cubit.cmd('block 2 attribute index 4 2670')   # rho 
cubit.cmd('block 2 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 
 

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH')  
 
# all files needed by SCOTCH are now in directory MESH 

