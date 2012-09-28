!=====================================================================
!
!               s p e c f e m 3 d  v e r s i o n  1 . 4
!               ---------------------------------------
!
!                 dimitri komatitsch and jeroen tromp
!    seismological laboratory - california institute of technology
!         (c) california institute of technology september 2006
!
! this program is free software; you can redistribute it and/or modify
! it under the terms of the gnu general public license as published by
! the free software foundation; either version 2 of the license, or
! (at your option) any later version.
!
! this program is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the
! gnu general public license for more details.
!
! you should have received a copy of the gnu general public license along
! with this program; if not, write to the free software foundation, inc.,
! 51 franklin street, fifth floor, boston, ma 02110-1301 usa.
!
!===============================================================================================================

! This module was written by:
! Percy Galvez , Jean-Paul Ampuero and Javier Ruiz
! based on fault_solver.f90

module fault_solver_kinematic

 implicit none  

 include 'constants.h'

 private

! outputs on selected fault nodes at every time step:
! slip, slip velocity, fault stresses
 type dataT_type
   integer                                    :: npoin
   integer, dimension(:), pointer             :: iglob
   real(kind=CUSTOM_REAL), dimension(:,:), pointer  :: d1,v1,t1,d2,v2,t2,t3
   character(len=70), dimension(:), pointer   :: name
 end type dataT_type

! DATAXZ_type used to read snapshots (temporal)
  type dataXZ_type
    real(kind=CUSTOM_REAL), dimension(:), pointer   :: d1, d2, v1, v2, & !Slip and Slip rate.
                                                       t1, t2, t3 !Tractions.
    real(kind=CUSTOM_REAL), dimension(:), pointer   :: xcoord,ycoord,zcoord  
    integer                                         :: npoin
  end type dataXZ_type

 type bc_kinflt_type
   private
   integer :: nspec,nglob
   real(kind=CUSTOM_REAL) :: dt
   real(kind=CUSTOM_REAL), dimension(:), pointer      :: B,invM1,invM2,Z
   real(kind=CUSTOM_REAL), dimension(:,:), pointer    :: T,D,V,coord
   real(kind=CUSTOM_REAL), dimension(:,:,:), pointer  :: R
   integer, dimension(:), pointer               :: ibulk1, ibulk2
   type(dataT_type)                             :: dataT
   type(dataXZ_type)                            :: dataXZ
   real(kind=CUSTOM_REAL) :: kin_dt
   integer  :: kin_it
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: v_kin_t1,v_kin_t2
 end type bc_kinflt_type

 type(bc_kinflt_type), allocatable, save        :: faults(:)

!Number of time steps defined by the user : NTOUT
 integer, save                :: NTOUT,NSNAP

 logical, save :: SIMULATION_TYPE_KIN = .false.
 
 public :: BC_KINFLT_init, BC_KINFLT_set_all, SIMULATION_TYPE_KIN


contains


!=====================================================================
! BC_KINFLT_init initializes kinematic faults 
!
! prname        fault database is read from file prname_fault_db.bin
! Minv          inverse mass matrix
! dt            global time step
!
subroutine BC_KINFLT_init(prname,Minv,DTglobal,nt)

 character(len=256), intent(in) :: prname ! 'proc***'
 real(kind=CUSTOM_REAL), intent(in) :: Minv(:)
 double precision, intent(in) :: DTglobal 
 integer, intent(in) :: nt

 real(kind=CUSTOM_REAL) :: dt
 integer :: iflt,ier,dummy_idfault
 integer :: nbfaults
  integer :: SIMULATION_TYPE
 character(len=256) :: filename
 integer, parameter :: IIN_PAR =151
 integer, parameter :: IIN_BIN =170
 real(kind=CUSTOM_REAL) :: DUMMY 

 NAMELIST / BEGIN_FAULT / dummy_idfault 

 dummy_idfault = 0

 open(unit=IIN_PAR,file='DATA/FAULT/Par_file_faults',status='old',iostat=ier)
 if( ier /= 0 ) then
   write(6,*) 'Have not found Par_file_faults: assume no faults' 
   return 
 endif

 dt = real(DTglobal)
 filename = prname(1:len_trim(prname))//'fault_db.bin'
 open(unit=IIN_BIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
 if( ier /= 0 ) stop 'Have not found proc*_fault_db.bin'
! WARNING TO DO: should be an MPI abort

  read(IIN_PAR,*) nbfaults
  if (nbfaults==0) return
  read(IIN_PAR,*)  ! eta
  read(IIN_PAR,*) SIMULATION_TYPE
  if ( SIMULATION_TYPE == 2 ) then
    SIMULATION_TYPE_KIN = .true.
    read(IIN_PAR,*) NTOUT
    read(IIN_PAR,*) NSNAP
    read(IIN_PAR,*) DUMMY
    read(IIN_PAR,*) DUMMY 
    read(IIN_BIN) nbfaults ! should be the same as in IIN_PAR
    allocate( faults(nbfaults) )
    do iflt=1,nbfaults
      read(IIN_PAR,nml=BEGIN_FAULT,end=100)
      call init_one_fault(faults(iflt),IIN_BIN,IIN_PAR,Minv,dt,nt,iflt)
    enddo 
  endif
  close(IIN_BIN)
  close(IIN_PAR)

  return
100 stop 'Did not find BEGIN_FAULT block #'
  ! WARNING TO DO: should be an MPI abort

end subroutine BC_KINFLT_init


!---------------------------------------------------------------------

subroutine init_one_fault(bc,IIN_BIN,IIN_PAR,Minv,dt,NT,iflt)

 type(bc_kinflt_type), intent(inout) :: bc
 real(kind=CUSTOM_REAL), intent(in)  :: Minv(:)
 integer, intent(in)                 :: IIN_BIN,IIN_PAR,NT,iflt
 real(kind=CUSTOM_REAL), intent(in)  :: dt

 real(kind=CUSTOM_REAL), dimension(:,:), allocatable   :: jacobian2Dw
 real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: normal
 integer, dimension(:,:), allocatable :: ibool1
 real(kind=CUSTOM_REAL) :: norm
 integer :: ij,k,e
 real(kind=CUSTOM_REAL), dimension(:), allocatable :: nx,ny,nz
 real(kind=CUSTOM_REAL) :: kindt

 NAMELIST / KINPAR / kindt

 read(IIN_BIN) bc%nspec,bc%nglob
 if (bc%nspec==0) return

 allocate( bc%ibulk1(bc%nglob) )
 allocate( bc%ibulk2(bc%nglob) )
 allocate( ibool1(NGLLSQUARE,bc%nspec) )
 allocate(normal(NDIM,NGLLSQUARE,bc%nspec))
 allocate(jacobian2Dw(NGLLSQUARE,bc%nspec))
 allocate(bc%coord(3,bc%nglob))

 read(IIN_BIN) ibool1
 read(IIN_BIN) jacobian2Dw
 read(IIN_BIN) normal
 read(IIN_BIN) bc%ibulk1
 read(IIN_BIN) bc%ibulk2
 read(IIN_BIN) bc%coord(1,:)
 read(IIN_BIN) bc%coord(2,:)
 read(IIN_BIN) bc%coord(3,:)
 bc%dt = dt

 allocate( bc%B(bc%nglob) ) 
 bc%B = 0e0_CUSTOM_REAL
 allocate( nx(bc%nglob),ny(bc%nglob),nz(bc%nglob) )
 nx = 0e0_CUSTOM_REAL
 ny = 0e0_CUSTOM_REAL
 nz = 0e0_CUSTOM_REAL
 do e=1,bc%nspec
   do ij = 1,NGLLSQUARE
     k = ibool1(ij,e)
     nx(k) = nx(k) + normal(1,ij,e)
     ny(k) = ny(k) + normal(2,ij,e)
     nz(k) = nz(k) + normal(3,ij,e)
     bc%B(k) = bc%B(k) + jacobian2Dw(ij,e)
   enddo
 enddo
 ! TO DO: assemble B and n across processors
 do k=1,bc%nglob
   norm = sqrt( nx(k)*nx(k) + ny(k)*ny(k) + nz(k)*nz(k) )
   nx(k) = nx(k) / norm
   ny(k) = ny(k) / norm 
   nz(k) = nz(k) / norm 
 enddo
 allocate( bc%R(3,3,bc%nglob) )
 call compute_R(bc%R,bc%nglob,nx,ny,nz)
 deallocate(nx,ny,nz)
! Needed in dA_Free = -K2*d2/M2 + K1*d1/M1
 allocate(bc%invM1(bc%nglob))
 allocate(bc%invM2(bc%nglob))
 bc%invM1 = Minv(bc%ibulk1)
 bc%invM2 = Minv(bc%ibulk2)

! Fault impedance, Z in :  Trac=T_Stick-Z*dV
!   Z = 1/( B1/M1 + B2/M2 ) / (0.5*dt)
! T_Stick = Z*Vfree traction as if the fault was stuck (no displ discontinuity) 
! NOTE: same Bi on both sides, see note above
 allocate(bc%Z(bc%nglob))
 bc%Z = 1.e0_CUSTOM_REAL/(0.5e0_CUSTOM_REAL*dt * bc%B *( bc%invM1 + bc%invM2 ))

 allocate(bc%T(3,bc%nglob))
 allocate(bc%D(3,bc%nglob))
 allocate(bc%V(3,bc%nglob))
 bc%T = 0e0_CUSTOM_REAL
 bc%D = 0e0_CUSTOM_REAL
 bc%V = 0e0_CUSTOM_REAL
! Dt between two loaded slip rates
 
 read(IIN_PAR,nml=KINPAR) 
 bc%kin_dt = kindt
 
 bc%kin_it=0
! Always have in memory the slip-rate model at two times, t1 and t2, 
! spatially interpolated in the spectral element grid
 allocate(bc%v_kin_t1(2,bc%nglob))
 allocate(bc%v_kin_t2(2,bc%nglob))
 bc%v_kin_t1 = 0e0_CUSTOM_REAL
 bc%v_kin_t2 = 0e0_CUSTOM_REAL

 call init_dataT(bc%dataT,bc%coord,bc%nglob,NT,iflt)
 call init_dataXZ(bc%dataXZ,bc%nglob)

end subroutine init_one_fault

!---------------------------------------------------------------------
subroutine compute_R(R,nglob,nx,ny,nz)

 integer :: nglob 
 real(kind=CUSTOM_REAL), intent(out) :: R(3,3,nglob)
 real(kind=CUSTOM_REAL), dimension(nglob), intent(in) :: nx,ny,nz

 real(kind=CUSTOM_REAL), dimension(nglob) :: sx,sy,sz,dx,dy,dz,norm

! Percy , defining fault directions (in concordance with SCEC conventions) . 
!   fault coordinates (s,d,n) = (1,2,3)
!   s = strike , d = dip , n = n. 
!   1 = strike , 2 = dip , 3 = n.  
    norm = sqrt(nx*nx+ny*ny)
    sx =  ny/norm  
    sy = -nx/norm     
    sz = 0.e0_CUSTOM_REAL  

    norm = sqrt(sy*sy*nz*nz+sx*sx*nz*nz+(sy*nx-ny*sx)*(nx*sy-ny*sx))
    dx = -sy*nz/norm
    dy =  sx*nz/norm
    dz = (sy*nx-ny*sx)/norm
!Percy, dz is always dipwards = -1/norm , because (nx*sy-ny*sx)= - 1 

    R(1,1,:)=sx
    R(1,2,:)=sy
    R(1,3,:)=sz
    R(2,1,:)=dx
    R(2,2,:)=dy
    R(2,3,:)=dz
    R(3,1,:)=nx
    R(3,2,:)=ny
    R(3,3,:)=nz
  

end subroutine compute_R


!=====================================================================
! adds boundary term Bt to Force array for each fault.
!
subroutine BC_KINFLT_set_all(F,Vel,Dis)

 real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Vel,Dis
 real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: F

 integer :: iflt

 if (.not. allocated(faults)) return
 do iflt=1,size(faults)
   if (faults(iflt)%nspec>0) call BC_KINFLT_set_single(faults(iflt),F,Vel,Dis,iflt)
 enddo 

end subroutine BC_KINFLT_set_all

!---------------------------------------------------------------------
subroutine BC_KINFLT_set_single(bc,MxA,V,D,iflt) 

 use specfem_par, only:it,NSTEP 

 real(kind=CUSTOM_REAL), intent(inout) :: MxA(:,:)
 type(bc_kinflt_type), intent(inout) :: bc
 real(kind=CUSTOM_REAL), intent(in) :: V(:,:),D(:,:)
 integer,intent(in) :: iflt
 integer :: it_kin,itime 
 real(kind=CUSTOM_REAL), dimension(3,bc%nglob) :: T
 real(kind=CUSTOM_REAL), dimension(3,bc%nglob) :: dD,dV,dA,dV_free
 real(kind=CUSTOM_REAL) :: t1,t2
 real(kind=CUSTOM_REAL) :: half_dt,time

 half_dt = 0.5e0_CUSTOM_REAL*bc%dt

! get predicted values
 dD = get_jump(bc,D) ! dD_predictor
 dV = get_jump(bc,V) ! dV_predictor
 dA = get_weighted_jump(bc,MxA) ! dA_free

! rotate to fault frame (tangent,normal)
! component 3 is normal to the fault
 dD = rotate(bc,dD,1)
 dV = rotate(bc,dV,1) 
 dA = rotate(bc,dA,1)   

! Time marching
 time = it*bc%dt
! Slip_rate step "it_kin"
 it_kin = bc%kin_it*nint(bc%kin_dt/bc%dt)
! (nint : fortran round (nearest whole number) , 
!  if nint(a)=0.5 then "a" get upper bound )

! Loading the next slipt_rate one ahead it.
! This is done in case bc%kin_dt 
! if (it_kin == it) it_kin=it_kin+1 ! 


!NOTE : it and it_kin is being used due to integers are exact numbers.
 if (it > it_kin) then

   print*, 'it :'
   print*, it
   print*, 'it_kin'
   print*, it_kin

   bc%kin_it = bc%kin_it +1
   bc%v_kin_t1 = bc%v_kin_t2
   print*, 'loading v_kin_t2'
   !Temporal : just for snapshots file names kin_dt=0.1 , dt=0.0001 
   !snapshot(100=itime).. : itime=kin_it*(kin_dt/dt)             
   itime = bc%kin_it*nint(bc%kin_dt/bc%dt)
   call load_vslip_snapshots(bc%dataXZ,itime,bc%nglob,iflt)
!   loading slip rates 
   bc%v_kin_t2(1,:)=bc%dataXZ%v1
   bc%v_kin_t2(2,:)=bc%dataXZ%v2
   
   !linear interpolation in time between t1 and t2
   !REMARK , bc%kin_dt is the delta "t" between two snapshots.
     
 endif

   t1 = (bc%kin_it-1) * bc%kin_dt
   t2 = bc%kin_it * bc%kin_dt

! Kinematic velocity_rate
! bc%V : Imposed apriori and read from slip rate snapshots (from time reversal)
!                linear interpolate between consecutive kinematic time steps.
!                V will be given each time step.
 bc%V(1,:) = ( (t2 - time)*bc%v_kin_t1(1,:) + (time - t1)*bc%v_kin_t2(1,:) )/ bc%kin_dt
 bc%V(2,:) = ( (t2 - time)*bc%v_kin_t1(2,:) + (time - t1)*bc%v_kin_t2(2,:) )/ bc%kin_dt

!dV_free = dV_predictor + (dt/2)*dA_free 
 dV_free(1,:) = dV(1,:)+half_dt*dA(1,:)
 dV_free(2,:) = dV(2,:)+half_dt*dA(2,:)
 dV_free(3,:) = dV(3,:)+half_dt*dA(3,:)

! T = Z*( dV_free - V) , V known apriori as input.
! CONVENTION : T(ibulk1)=T=-T(ibulk2)
 T(1,:) = bc%Z * ( dV_free(1,:) -bc%V(1,:) )
 T(2,:) = bc%Z * ( dV_free(2,:) -bc%V(2,:) )
 T(3,:) = bc%Z * ( dV_free(3,:) )

! Save tractions
 bc%T = T

! Update slip in fault frame
 bc%D = dD

! Rotate tractions back to (x,y,z) frame 
 T = rotate(bc,T,-1)

! Add boundary term B*T to M*a
 MxA(1,bc%ibulk1) = MxA(1,bc%ibulk1) + bc%B*T(1,:)
 MxA(2,bc%ibulk1) = MxA(2,bc%ibulk1) + bc%B*T(2,:)
 MxA(3,bc%ibulk1) = MxA(3,bc%ibulk1) + bc%B*T(3,:)

 MxA(1,bc%ibulk2) = MxA(1,bc%ibulk2) - bc%B*T(1,:)
 MxA(2,bc%ibulk2) = MxA(2,bc%ibulk2) - bc%B*T(2,:)
 MxA(3,bc%ibulk2) = MxA(3,bc%ibulk2) - bc%B*T(3,:)

!-- intermediate storage of outputs --
 call store_dataT(bc%dataT,bc%D,bc%V,bc%T,it)

!-- OUTPUTS --
! write dataT every NTOUT time steps or at the end of simulation
 if ( mod(it,NTOUT) == 0 .or. it==NSTEP) call SCEC_write_dataT(bc%dataT,bc%dt,it)
! write dataXZ every NSNAP time steps
! if ( mod(it,NSNAP) == 0) call write_dataXZ(bc,it,iflt)


end subroutine BC_KINFLT_set_single

!===============================================================
function get_jump(bc,v) result(dv)

 type(bc_kinflt_type), intent(in) :: bc
 real(kind=CUSTOM_REAL), intent(in) :: v(:,:)
 real(kind=CUSTOM_REAL) :: dv(3,bc%nglob)

! diference between side 2 and side 1 of fault nodes. dv
 dv(1,:) = v(1,bc%ibulk2)-v(1,bc%ibulk1)
 dv(2,:) = v(2,bc%ibulk2)-v(2,bc%ibulk1)
 dv(3,:) = v(3,bc%ibulk2)-v(3,bc%ibulk1)

end function get_jump

!---------------------------------------------------------------------
function get_weighted_jump(bc,f) result(da)

 type(bc_kinflt_type), intent(in) :: bc
 real(kind=CUSTOM_REAL), intent(in) :: f(:,:)
 real(kind=CUSTOM_REAL) :: da(3,bc%nglob)

! diference between side 2 and side 1 of fault nodes. M-1 * F
 da(1,:) = bc%invM2*f(1,bc%ibulk2)-bc%invM1*f(1,bc%ibulk1)
 da(2,:) = bc%invM2*f(2,bc%ibulk2)-bc%invM1*f(2,bc%ibulk1) 
 da(3,:) = bc%invM2*f(3,bc%ibulk2)-bc%invM1*f(3,bc%ibulk1)

end function get_weighted_jump

!----------------------------------------------------------------------
function rotate(bc,v,fb) result(vr)

 type(bc_kinflt_type), intent(in) :: bc
 real(kind=CUSTOM_REAL), intent(in) :: v(3,bc%nglob)
 integer, intent(in) :: fb
 real(kind=CUSTOM_REAL) :: vr(3,bc%nglob)

! Percy, tangential direction Vt, equation 7 of Pablo's notes in agreement with SPECFEM3D

! forward rotation
 if (fb==1) then
   vr(1,:) = v(1,:)*bc%R(1,1,:)+v(2,:)*bc%R(1,2,:)+v(3,:)*bc%R(1,3,:) ! vs
   vr(2,:) = v(1,:)*bc%R(2,1,:)+v(2,:)*bc%R(2,2,:)+v(3,:)*bc%R(2,3,:) ! vd
   vr(3,:) = v(1,:)*bc%R(3,1,:)+v(2,:)*bc%R(3,2,:)+v(3,:)*bc%R(3,3,:) ! vn

!  backward rotation
 else
   vr(1,:) = v(1,:)*bc%R(1,1,:)+v(2,:)*bc%R(2,1,:)+v(3,:)*bc%R(3,1,:)  !vx
   vr(2,:) = v(1,:)*bc%R(1,2,:)+v(2,:)*bc%R(2,2,:)+v(3,:)*bc%R(3,2,:)  !vy
   vr(3,:) = v(1,:)*bc%R(1,3,:)+v(2,:)*bc%R(2,3,:)+v(3,:)*bc%R(3,3,:)  !vz

 endif

end function rotate


!===============================================================
! OUTPUTS

subroutine init_dataT(DataT,coord,nglob,NT,iflt)
 ! NT = total number of time steps

 integer, intent(in) :: nglob,NT,iflt
 real(kind=CUSTOM_REAL), intent(in) :: coord(3,nglob)
 type (dataT_type), intent(out) :: DataT

 real(kind=CUSTOM_REAL) :: xtarget,ytarget,ztarget,dist,distkeep
 integer :: i, iglob , IIN, ier, jflt, np, k
 character(len=70) :: tmpname

!  1. read fault output coordinates from user file, 
!  2. define iglob: the fault global index of the node nearest to user
!     requested coordinate

 IIN = 251
 open(IIN,file='DATA/FAULT/FAULT_STATIONS',status='old',action='read',iostat=ier)
 read(IIN,*) np
 DataT%npoin =0
 do i=1,np
   read(IIN,*) xtarget,ytarget,ztarget,tmpname,jflt
   if (jflt==iflt) DataT%npoin = DataT%npoin +1
 enddo  
 close(IIN)

 if (DataT%npoin == 0) return

 allocate(DataT%iglob(DataT%npoin))
 allocate(DataT%name(DataT%npoin))

 open(IIN,file='DATA/FAULT/FAULT_STATIONS',status='old',action='read',iostat=ier)
 if( ier /= 0 ) stop 'error opening FAULT_STATIONS file'
 read(IIN,*) np
 k = 0
 do i=1,np
   read(IIN,*) xtarget,ytarget,ztarget,tmpname,jflt
   if (jflt/=iflt) cycle
   k = k+1
   DataT%name(k) = tmpname
  !search nearest node
   distkeep = huge(distkeep)

   do iglob=1,nglob
     dist = sqrt((coord(1,iglob)-xtarget)**2   &
          + (coord(2,iglob)-ytarget)**2 &
          + (coord(3,iglob)-ztarget)**2)  
     if (dist < distkeep) then
       distkeep = dist
       DataT%iglob(k) = iglob   
     endif 
   enddo
 enddo  

!  3. allocate arrays and set to zero
 allocate(DataT%d1(NT,DataT%npoin))
 allocate(DataT%v1(NT,DataT%npoin))
 allocate(DataT%t1(NT,DataT%npoin))
 allocate(DataT%d2(NT,DataT%npoin))
 allocate(DataT%v2(NT,DataT%npoin))
 allocate(DataT%t2(NT,DataT%npoin))
 allocate(DataT%t3(NT,DataT%npoin))
 DataT%d1 = 0e0_CUSTOM_REAL
 DataT%v1 = 0e0_CUSTOM_REAL
 DataT%t1 = 0e0_CUSTOM_REAL
 DataT%d2 = 0e0_CUSTOM_REAL
 DataT%v2 = 0e0_CUSTOM_REAL
 DataT%t2 = 0e0_CUSTOM_REAL
 DataT%t3 = 0e0_CUSTOM_REAL

 close(IIN)

end subroutine init_dataT


!---------------------------------------------------------------
subroutine init_dataXZ(dataXZ,nglob)

 type(dataXZ_type), intent(inout) :: dataXZ
 integer, intent(in) :: nglob

  allocate(dataXZ%v1(nglob))
  allocate(dataXZ%v2(nglob))
  allocate(dataXZ%xcoord(nglob))
  allocate(dataXZ%ycoord(nglob))
  allocate(dataXZ%zcoord(nglob))

  dataXZ%v1= 0e0_CUSTOM_REAL
  dataXZ%v2= 0e0_CUSTOM_REAL
  dataXZ%xcoord= 0e0_CUSTOM_REAL
  dataXZ%ycoord= 0e0_CUSTOM_REAL
  dataXZ%zcoord= 0e0_CUSTOM_REAL

end subroutine init_dataXZ


!---------------------------------------------------------------
subroutine store_dataT(dataT,d,v,t,itime)

 type(dataT_type), intent(inout) :: dataT
 real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: d,v,t
 integer, intent(in) :: itime

 integer :: i,k

 do i=1,dataT%npoin
   k = dataT%iglob(i)
   dataT%d1(itime,i) = d(1,k)
   dataT%d2(itime,i) = d(2,k)
   dataT%v1(itime,i) = v(1,k)
   dataT%v2(itime,i) = v(2,k)
   dataT%t1(itime,i) = t(1,k)
   dataT%t2(itime,i) = t(2,k)
   dataT%t3(itime,i) = t(3,k)
 enddo

end subroutine store_dataT


!-----------------------------------------------------------------

subroutine SCEC_write_dataT(dataT,DT,NT)

 type(dataT_type), intent(in) :: dataT
 real(kind=CUSTOM_REAL), intent(in) :: DT
 integer, intent(in) :: NT

 integer   :: i,k,IOUT

 IOUT = 121 !WARNING: not very robust. Could instead look for an available ID

do i=1,dataT%npoin

     open(IOUT,file='OUTPUT_FILES/'//trim(dataT%name(i))//'.dat',status='replace')
     write(IOUT,*) "% problem=TPV5"
     write(IOUT,*) "% author=Galvez, Ampuero, Nissen-Meyer"
     write(IOUT,*) "% date=2010/xx/xx"
     write(IOUT,*) "% code=SPECFEM3D_FAULT "
     write(IOUT,*) "% code_version=1.1"
     write(IOUT,*) "% element_size=100 m  (*4 GLL nodes)"
     write(IOUT,*) "% time_step=",DT
     write(IOUT,*) "% num_time_steps=",NT
     write(IOUT,*) "% location=",trim(dataT%name(i))
     write(IOUT,*) "% Time series in 8 column of E15.7"
     write(IOUT,*) "% Column #1 = Time (s)"
     write(IOUT,*) "% Column #2 = horizontal right-lateral slip (m)"
     write(IOUT,*) "% Column #3 = horizontal right-lateral slip rate (m/s)"
     write(IOUT,*) "% Column #4 = horizontal right-lateral shear stress (MPa)"
     write(IOUT,*) "% Column #5 = vertical up-dip slip (m)"
     write(IOUT,*) "% Column #6 = vertical up-dip slip rate (m/s)"
     write(IOUT,*) "% Column #7 = vertical up-dip shear stress (MPa)"
     write(IOUT,*) "% Column #8 = normal stress (MPa)"
     write(IOUT,*) "%"
     write(IOUT,*) "% The line below lists the names of the data fields:"
     write(IOUT,*) "%t  h-slip  h-slip-rate  h-shear-stress v-slip v-slip-rate v-shear-stress n-stress"
     write(IOUT,*) "%"
     write(IOUT,*) "% Here is the time-series data."
     do k=1,NT
       write(IOUT,'(8(E15.7))') k*DT, dataT%d1(k,i), dataT%v1(k,i), dataT%t1(k,i)/1.0e6_CUSTOM_REAL, &
                                        dataT%d2(k,i), dataT%v2(k,i), dataT%t2(k,i)/1.0e6_CUSTOM_REAL, &
                                        dataT%t3(k,i)/1.0e6_CUSTOM_REAL
     enddo
     close(IOUT)
 enddo

end subroutine SCEC_write_dataT


!---------------------------------------------------------------
!LOAD_VSLIP_SNAPSHOTS(v,dataXZ,itime,coord,npoin,nglob,iflt)  
!Loading slip velocity from snapshots.
!   INPUT  itime : iteration time
!          coord : Receivers coordinates
!          npoin : number of Receivers.
!          nglob : number of gll points along the fault.
!          dataXZ : slip rate .
!          iflt : number of faults.

!   OUTPUT v : slip rate on receivers.
 
subroutine load_vslip_snapshots(dataXZ,itime,nglob,iflt)  

  integer, intent(in) :: itime,nglob,iflt
  type(dataXZ_type), intent(inout) :: dataXZ
  character(len=70) :: filename
  integer :: IIN_BIN,ier,IOUT

  IIN_BIN=101
  IOUT = 102

  write(filename,"('OUTPUT_FILES/Snapshot',I0,'_F',I0,'.bin')") itime,iflt
  print*, trim(filename)

  open(unit=IIN_BIN, file= trim(filename), status='old', form='formatted',&
        action='read',iostat=ier)
!  COMPILLERS WRITE BINARY OUTPUTS IN DIFFERENT FORMATS !!!!!!!!!! 
!  open(unit=IIN_BIN, file= trim(filename), status='old', form='unformatted',&
!        action='read',iostat=ier)
!  if( ier /= 0 ) stop 'Snapshots have been found'
 
  read(IIN_BIN,"(5F24.15)") dataXZ%xcoord,dataXZ%ycoord,dataXZ%zcoord,dataXZ%v1,dataXZ%v2

!  read(IOUT) dataXZ%xcoord
!  read(IOUT) dataXZ%ycoord
!  read(IOUT) dataXZ%zcoord
!  write(IOUT) dataXZ%d1
!  write(IOUT) dataXZ%d2
!  read(IOUT) dataXZ%v1
!  read(IOUT) dataXZ%v2
!  write(IOUT) dataXZ%t1
!  write(IOUT) dataXZ%t2
!  write(IOUT) dataXZ%t3
!  write(IOUT) dataXZ%sta
!  write(IOUT) dataXZ%stg
!  write(IOUT) dataXZ%tRUP
!  write(IOUT) dataXZ%tPZ
  close(IOUT)

  close(IIN_BIN)

end subroutine load_vslip_snapshots
!---------------------------------------------------------------

end module fault_solver_kinematic

