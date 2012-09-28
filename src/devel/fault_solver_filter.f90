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
! Percy Galvez, Jean-Paul Ampuero and Tarje Nissen-Meyer

module fault_solver

  implicit none  

  include 'constants.h'

  private

 ! outputs on selected fault nodes at every time step:
 ! slip, slip velocity, fault stresses
  type dataT_type
    integer                                    :: npoin
    integer, dimension(:), pointer             :: iglob   ! on-fault global index of output nodes
    real(kind=CUSTOM_REAL), dimension(:,:), pointer  :: d1,v1,t1,d2,v2,t2,t3
    character(len=70), dimension(:), pointer   :: name
  end type dataT_type

  
 ! outputs at selected times for all fault nodes:
 ! strength, state, slip, slip velocity, fault stresses, rupture time, process zone time
 ! rupture time = first time when slip velocity = threshold V_RUPT (defined below)
 ! process zone time = first time when slip = Dc
  type dataXZ_type
    real(kind=CUSTOM_REAL), dimension(:), pointer   :: stg, sta, d1, d2, v1, v2, & 
                                                       t1, t2, t3, tRUP,tPZ
    real(kind=CUSTOM_REAL), dimension(:), pointer   :: xcoord,ycoord,zcoord  
    integer                                         :: npoin
  end type dataXZ_type

  type swf_type
    private
    integer :: kind
    logical :: healing = .false.
    real(kind=CUSTOM_REAL), dimension(:), pointer   :: Dc=>null(), Dr=>null(), mus=>null(), &
                                                       mud=>null(),theta=>null(),idrop=>null()
                                                       ! idrop for multiple drops in friction.
  end type swf_type

  type bc_dynflt_type
    private
    integer :: nspec,nglob
    real(kind=CUSTOM_REAL), dimension(:,:), pointer    :: T0,T,V,D
    real(kind=CUSTOM_REAL), dimension(:,:), pointer    :: coord 
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer  :: R
    real(kind=CUSTOM_REAL), dimension(:), pointer      :: MU,B,invM1,invM2,Z  
    real(kind=CUSTOM_REAL) :: dt
    integer, dimension(:), pointer               :: ibulk1, ibulk2  
    type(swf_type), pointer                      :: swf => null(), swfi => null()   ! swfi for multiple friction drops
    logical                                      :: allow_opening = .false.         ! default : do not allow opening
    type(dataT_type)                             :: dataT
    type(dataXZ_type)                            :: dataXZ
  end type bc_dynflt_type

  type(bc_dynflt_type), allocatable, save        :: faults(:)

 !slip velocity threshold for healing
 !WARNING: not very robust
  real(kind=CUSTOM_REAL), save       :: V_HEALING 

 !slip velocity threshold for definition of rupture front
  real(kind=CUSTOM_REAL), save       :: V_RUPT 

 !Number of time steps defined by the user : NTOUT
  integer, save                :: NTOUT,NSNAP
 !Type of input for initial fault Parameters 
  integer, save                :: INPUT_FAULT_PAR_TYPE,NDROPS ! NDROPS : number of drops in linear  
                                                              !          weak slip friction coeficient 
  logical, save :: SIMULATION_TYPE_DYN = .false.

  real(kind=CUSTOM_REAL), allocatable, save :: Kelvin_Voigt_eta(:)
  
  public :: BC_DYNFLT_init, BC_DYNFLT_set3d_all, Kelvin_Voigt_eta, &
            SIMULATION_TYPE_DYN


contains


!=====================================================================
! BC_DYNFLT_init initializes dynamic faults 
!
! prname        fault database is read from file prname_fault_db.bin
! Minv          inverse mass matrix
! dt            global time step
!
  subroutine BC_DYNFLT_init(prname,Minv,DTglobal,nt)

  character(len=256), intent(in) :: prname ! 'proc***'
  real(kind=CUSTOM_REAL), intent(in) :: Minv(:)
  double precision, intent(in) :: DTglobal 
  integer, intent(in) :: nt

  real(kind=CUSTOM_REAL) :: dt
  integer :: iflt,ier,dummy_idfault
  integer :: nbfaults
  integer :: size_Kelvin_Voigt
  integer :: SIMULATION_TYPE
  character(len=256) :: filename
  integer, parameter :: IIN_PAR =151
  integer, parameter :: IIN_BIN =170
  integer, parameter :: IIN_FAL =201

  NAMELIST / BEGIN_FAULT / dummy_idfault 

  dummy_idfault = 0

  open(unit=IIN_PAR,file='DATA/FAULT/Par_file_faults',status='old',iostat=ier)
  if( ier /= 0 ) then
    write(6,*) 'File Par_file_faults not found: assume no faults'
    close(IIN_PAR) 
    return 
  endif

  dt = real(DTglobal)
  filename = prname(1:len_trim(prname))//'fault_db.bin'
  open(unit=IIN_BIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    write(6,*) 'File ',trim(filename),' not found. Abort' 
    stop
  endif
! WARNING TO DO: should be an MPI abort
  
  read(IIN_PAR,*) nbfaults
  if (nbfaults==0) return
! Reading etas of each fault
  do iflt = 1,nbfaults
    read(IIN_PAR,*) ! etas
  enddo
  read(IIN_PAR,*) SIMULATION_TYPE
  if ( SIMULATION_TYPE /= 1 ) then
    close(IIN_BIN)
    close(IIN_PAR)
    return 
  endif
  SIMULATION_TYPE_DYN = .true.
  read(IIN_PAR,*) NTOUT
  read(IIN_PAR,*) NSNAP 
  read(IIN_PAR,*) V_HEALING
  read(IIN_PAR,*) V_RUPT
  read(IIN_PAR,*) INPUT_FAULT_PAR_TYPE
  read(IIN_PAR,*) NDROPS 
 
  read(IIN_BIN) nbfaults ! should be the same as in IIN_PAR
  allocate( faults(nbfaults) )
  do iflt=1,nbfaults
    read(IIN_PAR,nml=BEGIN_FAULT,end=100)
    call init_one_fault(faults(iflt),IIN_BIN,IIN_PAR,IIN_FAL,Minv,dt,nt,iflt)
  enddo 
  close(IIN_BIN)
  close(IIN_PAR)
  
  filename = prname(1:len_trim(prname))//'Kelvin_voigt_eta.bin'
  open(unit=IIN_BIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    write(6,*) 'File ',trim(filename),' not found. Abort' 
    stop
  endif
  read(IIN_BIN) size_Kelvin_Voigt
  if (size_Kelvin_Voigt > 0) then
    allocate(Kelvin_Voigt_eta(size_Kelvin_Voigt))
    read(IIN_BIN) Kelvin_Voigt_eta
  endif
  close(IIN_BIN)
  return
100 stop 'Did not find BEGIN_FAULT block #'
   ! WARNING TO DO: should be an MPI abort

end subroutine BC_DYNFLT_init


!---------------------------------------------------------------------

  subroutine init_one_fault(bc,IIN_BIN,IIN_PAR,IIN_FAL,Minv,dt,NT,iflt)
  
  type(bc_dynflt_type), intent(inout) :: bc
  real(kind=CUSTOM_REAL), intent(in)  :: Minv(:)
  integer, intent(in)                 :: IIN_BIN,IIN_PAR,IIN_FAL,NT,iflt

  real(kind=CUSTOM_REAL), intent(in)  :: dt

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable   :: jacobian2Dw
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: normal
  integer, dimension(:,:), allocatable :: ibool1
  real(kind=CUSTOM_REAL) :: norm
  real(kind=CUSTOM_REAL) :: S1,S2,S3
  integer :: n1,n2,n3
  real(kind=CUSTOM_REAL) :: mus,mud,dc,dr
  integer :: nmus,nmud,ndc,ndr,ij,k,e
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: nx,ny,nz


  NAMELIST / INIT_STRESS / S1,S2,S3,n1,n2,n3
  NAMELIST / SWF / mus,mud,dc,dr,nmus,nmud,ndc,ndr

  read(IIN_BIN) bc%nspec,bc%nglob
  if (bc%nspec==0) return
  allocate( bc%ibulk1(bc%nglob) )
  allocate( bc%ibulk2(bc%nglob) )
  allocate( ibool1(NGLLSQUARE,bc%nspec) )
  allocate(normal(NDIM,NGLLSQUARE,bc%nspec))
  allocate(jacobian2Dw(NGLLSQUARE,bc%nspec))
  
  allocate(bc%coord(3,(bc%nglob)))
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
  do k=1,bc%nglob
    norm = sqrt( nx(k)*nx(k) + ny(k)*ny(k) + nz(k)*nz(k) )
    nx(k) = nx(k) / norm
    ny(k) = ny(k) / norm 
    nz(k) = nz(k) / norm 
  enddo

  allocate( bc%R(3,3,bc%nglob) )
  call compute_R(bc%R,bc%nglob,nx,ny,nz)

! Needed in dA_Free = -K2*d2/M2 + K1*d1/M1
  allocate(bc%invM1(bc%nglob))
  allocate(bc%invM2(bc%nglob))
  bc%invM1 = Minv(bc%ibulk1)
  bc%invM2 = Minv(bc%ibulk2)

! Fault impedance, Z in :  Trac=T_Stick-Z*dV
!   Z = 1/( B1/M1 + B2/M2 ) / (0.5*dt)
! T_stick = Z*Vfree traction as if the fault was stuck (no displ discontinuity) 
! NOTE: same Bi on both sides, see note above
  allocate(bc%Z(bc%nglob))
  bc%Z = 1.e0_CUSTOM_REAL/(0.5e0_CUSTOM_REAL*dt * bc%B *( bc%invM1 + bc%invM2 ))

  allocate(bc%T(3,bc%nglob))
  allocate(bc%D(3,bc%nglob))
  allocate(bc%V(3,bc%nglob))
  bc%T = 0e0_CUSTOM_REAL
  bc%D = 0e0_CUSTOM_REAL
  bc%V = 0e0_CUSTOM_REAL

! Set initial fault stresses
  allocate(bc%T0(3,bc%nglob))
  S1 = 0e0_CUSTOM_REAL
  S2 = 0e0_CUSTOM_REAL
  S3 = 0e0_CUSTOM_REAL
  n1=0
  n2=0
  n3=0
 
! WARNING : Quick and dirty free surface condition at z=0 
!  do k=1,bc%nglob  
!    if (abs(bc%zcoord(k)-0.e0_CUSTOM_REAL) <= SMALLVAL) bc%T0(2,k) = 0
!  end do 
! Set friction parameters and initialize friction variables
  allocate( bc%swf )
  allocate( bc%swf%mus(bc%nglob) )
  allocate( bc%swf%mud(bc%nglob) )
  allocate( bc%swf%Dc(bc%nglob) )
  allocate( bc%swf%Dr(bc%nglob) )
  allocate( bc%swf%theta(bc%nglob) )
! Set friction parameters and initialize friction variables updated
  allocate( bc%swfi )
  allocate( bc%swfi%mus(bc%nglob) )
  allocate( bc%swfi%mud(bc%nglob) )
  allocate( bc%swfi%theta(bc%nglob) )
  allocate( bc%swfi%idrop(bc%nglob) )
! WARNING: if V_HEALING is negative we turn off healing
  bc%swf%healing = (V_HEALING > 0e0_CUSTOM_REAL)

  mus = 0.6e0_CUSTOM_REAL 
  mud = 0.1e0_CUSTOM_REAL 
  dc = 1e0_CUSTOM_REAL
  dr = 10000e0_CUSTOM_REAL
  nmus = 0
  nmud = 0
  ndc  = 0
  ndr  = 0

  bc%swf%theta  = 0e0_CUSTOM_REAL
  bc%swfi%theta = 0e0_CUSTOM_REAL
  bc%swfi%idrop = 1

  allocate(bc%MU(bc%nglob))
 
  Select case (INPUT_FAULT_PAR_TYPE)
    case (1) 
      call init_2d_distribution_random(bc%coord,bc%T0,bc%swf%mus,bc%swf%mud, &
                                       bc%swf%Dc,bc%swf%Dr,IIN_FAL,iflt) 
    case (2)
      read(IIN_PAR, nml=INIT_STRESS)
      bc%T0(1,:) = S1
      bc%T0(2,:) = S2
      bc%T0(3,:) = S3
      call init_2d_distribution(bc%T0(1,:),bc%coord,IIN_PAR,n1) 
      call init_2d_distribution(bc%T0(2,:),bc%coord,IIN_PAR,n2) 
      call init_2d_distribution(bc%T0(3,:),bc%coord,IIN_PAR,n3) 
      read(IIN_PAR, nml=SWF)
      bc%swf%mus = mus
      bc%swf%mud = mud
      bc%swf%Dc  = dc
      bc%swf%Dr  = dr
      call init_2d_distribution(bc%swf%mus,bc%coord,IIN_PAR,nmus)
      call init_2d_distribution(bc%swf%mud,bc%coord,IIN_PAR,nmud) 
      call init_2d_distribution(bc%swf%Dc ,bc%coord,IIN_PAR,ndc)
      call init_2d_distribution(bc%swf%Dr ,bc%coord,IIN_PAR,ndr)
  End select
 
  bc%swfi%mus = bc%swf%mus
  bc%swfi%mud = bc%swf%mud 

  bc%MU = swf_mu(bc%swf)
  call init_dataT(bc%dataT,bc%coord,bc%nglob,NT,iflt)
  call init_dataXZ(bc%dataXZ,bc,bc%nglob)

!----  Initial Fault Parameters ------------- 
  call write_init_fault_par(bc,iflt)

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

!---------------------------------------------------------------------
! Adds a value to a fault parameter inside loaded from a fault paramters file
! Input x,y,z T01 T02 T03 mus mud Dc

  subroutine init_2d_distribution_random(coord,T0,mus,mud,Dc,Dr,IIN_FAUL,iflt)
  real(kind=CUSTOM_REAL), intent(inout) :: T0(:,:)
  real(kind=CUSTOM_REAL), intent(inout) :: mus(:),mud(:),Dc(:),Dr(:)
  integer, intent(in) :: IIN_FAUL,iflt
  real(kind=CUSTOM_REAL), intent(in) :: coord(:,:)
  
  character(len=256) :: filepar 
  integer            :: fault_nodes,ier
!  integer            :: inode
!  real(kind=CUSTOM_REAL)   :: Ti0(3,:)
!  real(kind=CUSTOM_REAL)   :: imus(:),imud(:),iDc(:),iDr(:)
!  integer                  :: icoord(3,:)
  
  write(filepar,"('DATA/FAULT/Par_init_faults_F',I0)") iflt
  
  open(unit=IIN_FAUL,file=trim(filepar),status='old',action='read',iostat=ier)
  if (ier==0) then 
    write(6,*) 'Initial fault parameters will be taken from Par_init_faults_par'     
  else
    write(6,*) 'Par_init_fault_par does not exist, '
    write(6,*) 'Initial fault parameters missing'
    stop  
  end if
  
  read(IIN_FAUL,*) fault_nodes
  if (fault_nodes /= size(coord,2)) then 
      write(6,*) 'number of faults nodes do not coincide '
      stop
  end if  
! WARNING
! Considering same index ordering than x,y,z,t0,mus,mud,dc,dr , however 
! this is a particular case , it needs to be changed
  read(IIN_FAUL,*)  T0(1,:)
  read(IIN_FAUL,*)  T0(2,:)
  read(IIN_FAUL,*)  T0(3,:)
  read(IIN_FAUL,*)  mus(:)
  read(IIN_FAUL,*)  mud(:)
  read(IIN_FAUL,*)  Dc(:)
  read(IIN_FAUL,*)  Dr(:) 

!  read(IIN_FAUL,*)  icoord(1,:),icoord(2,:),icoord(3,:),Ti0(1,:),Ti0(2,:),Ti0(3,:),&
!                    imus(:),imud(:),iDc(:),iDr(:) 
!  do inode = 1,fault_nodes
!    where ( (coord(1,:)==icoord(1,inode)).and.&
!            (coord(2,:)==icoord(2,inode)).and.&
!            (coord(3,:)==icoord(3,inode)) )
!             T0  = Ti0
!             mus = imus
!             mud = imud
!             Dc  = iDc
!             Dr  = iDr 
!    endwhere  
!  enddo
  close(IIN_FAUL)
  end subroutine init_2d_distribution_random


!---------------------------------------------------------------------
! adds a value to a fault parameter inside an area with prescribed shape
  subroutine init_2d_distribution(a,coord,iin,n)

  real(kind=CUSTOM_REAL), intent(inout) :: a(:)
  real(kind=CUSTOM_REAL), intent(in) :: coord(:,:)
  integer, intent(in) :: iin,n

  real(kind=CUSTOM_REAL) :: b(size(a))
  character(len=20) :: shape
  real(kind=CUSTOM_REAL) :: val,valh, xc, yc, zc, r, l, lx,ly,lz
  integer :: i

  NAMELIST / DIST2D / shape, val,valh, xc, yc, zc, r, l, lx,ly,lz

  if (n==0) return   
  
  do i=1,n
    shape = ''
    xc = 0e0_CUSTOM_REAL
    yc = 0e0_CUSTOM_REAL
    zc = 0e0_CUSTOM_REAL
    r = 0e0_CUSTOM_REAL
    l = 0e0_CUSTOM_REAL
    lx = 0e0_CUSTOM_REAL
    ly = 0e0_CUSTOM_REAL
    lz = 0e0_CUSTOM_REAL
    valh  = 0e0_CUSTOM_REAL

    read(iin,DIST2D)
    select case(shape)
      case ('circle')
        b = heaviside(r-sqrt((coord(1,:)-xc)**2 + &
                             (coord(2,:)-yc)**2 + (coord(3,:)-zc)**2)) * val
      case ('ellipse')
        b = heaviside(1e0_CUSTOM_REAL-sqrt((coord(1,:)-xc)**2/lx**2 + &
                                           (coord(2,:)-yc)**2/ly**2 + (coord(3,:)-zc)**2/lz**2)) * val
      case ('square')
        b = heaviside((l/2._CUSTOM_REAL)-abs(coord(1,:)-xc)+SMALLVAL) * & 
            heaviside((l/2._CUSTOM_REAL)-abs(coord(2,:)-yc)+SMALLVAL) * & 
            heaviside((l/2._CUSTOM_REAL)-abs(coord(3,:)-zc)+SMALLVAL) * &
            val
      case ('cilinder')
        b = heaviside(r - sqrt((coord(1,:)-xc)**2 + (coord(2,:)-yc)**2)) * &
            heaviside((lz/2._CUSTOM_REAL)-abs(coord(3,:)-zc)+SMALLVAL)   * & 
            val
      case ('cilinder-taper')
        b = heaviside(r - sqrt((coord(1,:)-xc)**2 + (coord(2,:)-yc)**2)) * &
            heaviside((lz/2._CUSTOM_REAL)-abs(coord(3,:)-zc)+SMALLVAL)   * & 
            (val + ( coord(3,:) - zc + lz/2._CUSTOM_REAL ) * ((valh-val)/lz))
      case ('rectangle')
        b = heaviside((lx/2._CUSTOM_REAL)-abs(coord(1,:)-xc)+SMALLVAL) * &
            heaviside((ly/2._CUSTOM_REAL)-abs(coord(2,:)-yc)+SMALLVAL) * &
            heaviside((lz/2._CUSTOM_REAL)-abs(coord(3,:)-zc)+SMALLVAL) * &
            val
      case ('rectangle-taper')
        b = heaviside((lx/2._CUSTOM_REAL)-abs(coord(1,:)-xc)+SMALLVAL) * &
            heaviside((ly/2._CUSTOM_REAL)-abs(coord(2,:)-yc)+SMALLVAL) * &
            heaviside((lz/2._CUSTOM_REAL)-abs(coord(3,:)-zc)+SMALLVAL) * &
            (val + ( coord(3,:) - zc + lz/2._CUSTOM_REAL ) * ((valh-val)/lz))
      case default
        stop 'bc_dynflt_3d::init_2d_distribution:: unknown shape'
    end select
   
    where (b /= 0e0_CUSTOM_REAL) a = b
  enddo
    
  end subroutine init_2d_distribution

!---------------------------------------------------------------------
  elemental function heaviside(x)

  real(kind=CUSTOM_REAL), intent(in) :: x
  real(kind=CUSTOM_REAL) :: heaviside

  if (x>=0e0_CUSTOM_REAL) then
    heaviside = 1e0_CUSTOM_REAL
  else
    heaviside = 0e0_CUSTOM_REAL
  endif

  end function heaviside

!=====================================================================
! adds boundary term Bt into Force array for each fault.
!
  subroutine bc_dynflt_set3d_all(F,Vel,Dis)

  real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Vel,Dis
  real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: F

  integer :: iflt

  if (.not. allocated(faults)) return
  do iflt=1,size(faults)
    if (faults(iflt)%nspec>0) call BC_DYNFLT_set3d(faults(iflt),F,Vel,Dis,iflt)
  enddo 
   
  end subroutine bc_dynflt_set3d_all

!---------------------------------------------------------------------
  subroutine BC_DYNFLT_set3d(bc,MxA,V,D,iflt) 
  
  use specfem_par, only:it,NSTEP 

  real(kind=CUSTOM_REAL), intent(inout) :: MxA(:,:)
  type(bc_dynflt_type), intent(inout) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: V(:,:),D(:,:)
  integer,intent(in) :: iflt

  real(kind=CUSTOM_REAL), dimension(bc%nglob) :: strength
  real(kind=CUSTOM_REAL), dimension(3,bc%nglob) :: T
  real(kind=CUSTOM_REAL), dimension(bc%nglob) :: t1,t2,tnorm,tnew
  real(kind=CUSTOM_REAL), dimension(3,bc%nglob) :: dD,dV,dA
  real(kind=CUSTOM_REAL), dimension(bc%nglob) :: theta_old, Vnorm, Vnorm_old
  real(kind=CUSTOM_REAL) :: half_dt
!  integer :: k 
 
  half_dt = 0.5e0_CUSTOM_REAL*bc%dt
  theta_old = bc%swf%theta
  Vnorm_old = sqrt(bc%V(1,:)*bc%V(1,:)+bc%V(2,:)*bc%V(2,:))

! get predicted values
  dD = get_jump(bc,D) ! dD_predictor
  dV = get_jump(bc,V) ! dV_predictor
  dA = get_weighted_jump(bc,MxA) ! dA_free

! rotate to fault frame (tangent,normal)
! component 3 is normal to the fault
  dD = rotate(bc,dD,1)
  dV = rotate(bc,dV,1) 
  dA = rotate(bc,dA,1)   

! T_stick
 T(1,:) = bc%Z * ( dV(1,:) + half_dt*dA(1,:) )
 T(2,:) = bc%Z * ( dV(2,:) + half_dt*dA(2,:) )
 T(3,:) = bc%Z * ( dV(3,:) + half_dt*dA(3,:) )

!Warning : dirty particular free surface condition z = 0. 
!  where (bc%zcoord(:) > - SMALLVAL) T(2,:) = 0
! do k=1,bc%nglob  
!   if (abs(bc%zcoord(k)-0.e0_CUSTOM_REAL) < SMALLVAL) T(2,k) = 0.e0_CUSTOM_REAL
! end do 

! add initial stress
  T = T + bc%T0
 
! Solve for normal stress (negative is compressive)
  ! Opening implies free stress
   if (bc%allow_opening) T(3,:) = min(T(3,:),0.e0_CUSTOM_REAL) 

! Update slip weakening friction:
! Update slip state variable
! WARNING: during opening the friction state variable should not evolve
  call swf_update_state(bc%D,dD,bc%V,bc%swf)
 
! Update friction coeficient
! bc%MU = swf_mu(bc%swf) 

! Update Multi_step slip weakening 
  bc%MU = swf_multi_mu(bc%swf,bc%swfi,bc%T0,theta_old) 

! combined with time-weakening for nucleation
!  if (associated(bc%twf)) bc%MU = min( bc%MU, twf_mu(bc%twf,bc%coord,time) )

! Update strength
  strength = -bc%MU * min(T(3,:),0.e0_CUSTOM_REAL)

! Solve for shear stress
  tnorm = sqrt( T(1,:)*T(1,:) + T(2,:)*T(2,:))
  tnorm = max(tnorm,1e0_CUSTOM_REAL)
  t1 = T(1,:)/tnorm
  t2 = T(2,:)/tnorm

  tnew = min(tnorm,strength) 
  T(1,:) = tnew * t1
  T(2,:) = tnew * t2

! Save total tractions
  bc%T = T

! Subtract initial stress
  T = T - bc%T0

! Update slip acceleration da=da_free-T/(0.5*dt*Z)
  dA(1,:) = dA(1,:) - T(1,:)/(bc%Z*half_dt)
  dA(2,:) = dA(2,:) - T(2,:)/(bc%Z*half_dt)
  dA(3,:) = dA(3,:) - T(3,:)/(bc%Z*half_dt)
   
! Update slip and slip rate, in fault frame
  bc%D = dD
  bc%V = dV + half_dt*dA

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
  Vnorm = sqrt(bc%V(1,:)*bc%V(1,:)+bc%V(2,:)*bc%V(2,:))
  call store_dataXZ(bc%dataXZ, strength, theta_old, bc%swf%theta, bc%swf%dc, &
                    Vnorm_old, Vnorm, it*bc%dt,bc%dt)
  call store_dataT(bc%dataT,bc%D,bc%V,bc%T,it)


!-- outputs --
! write dataT every NTOUT time step or at the end of simulation
  if ( mod(it,NTOUT) == 0 .or. it==NSTEP) call SCEC_write_dataT(bc%dataT,bc%dt,it)
! write dataXZ every NSNAP time step
  if ( mod(it,NSNAP) == 0) call write_dataXZ(bc%dataXZ,it,iflt)
  if ( it == NSTEP) call SCEC_Write_RuptureTime(bc%dataXZ,bc%dt,NSTEP,iflt)

  end subroutine BC_DYNFLT_set3d

!===============================================================
 function get_jump (bc,v) result(dv)

  type(bc_dynflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: v(:,:)
  real(kind=CUSTOM_REAL) :: dv(3,bc%nglob)

! diference between side 2 and side 1 of fault nodes. dv
    dv(1,:) = v(1,bc%ibulk2)-v(1,bc%ibulk1)
    dv(2,:) = v(2,bc%ibulk2)-v(2,bc%ibulk1)
    dv(3,:) = v(3,bc%ibulk2)-v(3,bc%ibulk1)
    
  end function get_jump

!---------------------------------------------------------------------
  function get_weighted_jump (bc,f) result(da)

    type(bc_dynflt_type), intent(in) :: bc
    real(kind=CUSTOM_REAL), intent(in) :: f(:,:)

    real(kind=CUSTOM_REAL) :: da(3,bc%nglob)

! diference between side 2 and side 1 of fault nodes. M-1 * F
     da(1,:) = bc%invM2*f(1,bc%ibulk2)-bc%invM1*f(1,bc%ibulk1)
     da(2,:) = bc%invM2*f(2,bc%ibulk2)-bc%invM1*f(2,bc%ibulk1) 
     da(3,:) = bc%invM2*f(3,bc%ibulk2)-bc%invM1*f(3,bc%ibulk1)
  
  end function get_weighted_jump

!----------------------------------------------------------------------
  function rotate(bc,v,fb) result(vr)

  type(bc_dynflt_type), intent(in) :: bc
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


!=====================================================================
  subroutine swf_update_state(dold,dnew,vold,f)

  real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: vold,dold,dnew
  type(swf_type), intent(inout) :: f

  real(kind=CUSTOM_REAL) :: vnorm
  integer :: k,npoin

  f%theta = f%theta + sqrt( (dold(1,:)-dnew(1,:))**2 + (dold(2,:)-dnew(2,:))**2 )

  if (f%healing) then
    npoin = size(vold,2) 
    do k=1,npoin
      vnorm = sqrt(vold(1,k)**2 + vold(2,k)**2)
      if (vnorm<V_HEALING) f%theta(k) = 0e0_CUSTOM_REAL
    enddo
  endif
  end subroutine swf_update_state


!=====================================================================
! Friction coefficient
  function swf_mu(f) result(mu)

  type(swf_type), intent(in) :: f
  real(kind=CUSTOM_REAL) :: mu(size(f%theta))

 !-- linear slip weakening:

    mu = f%mus -(f%mus-f%mud)/f%dc *f%theta
    mu = max( mu, f%mud)
 
  end function swf_mu


!=====================================================================
! Stress Drop 

  function stress_drop(T0,mud) result(str_drop)
   
   real(kind=custom_real),dimension(:),intent(in)    :: mud
   real(kind=custom_real),dimension(3,size(mud)),intent(in)  :: T0

   real(kind=custom_real),dimension(size(mud))        :: str_drop,tnorm
   

    tnorm  = sqrt(T0(1,:)**2 + T0(2,:)**2)
!    str_drop = tnorm - mud*abs(T0(3,:))-c
    str_drop = tnorm - mud*abs(T0(3,:))

  end function stress_drop 


!=====================================================================
! Friction coefficient

  function swf_multi_mu(f,fi,T0,theta_old) result(mu)

  type(swf_type), intent(in)                                    :: f
  real(kind=custom_real),dimension(3,size(f%theta)),intent(in)  :: T0
  real(kind=custom_real),dimension(size(f%theta)),intent(in)    :: theta_old

  real(kind=custom_real),dimension(size(f%theta))               :: mu,u_drop

  type(swf_type), intent(inout)                                 :: fi


!-- multi linear slip weakening
!   theta : slip

   u_drop = stress_drop(T0,f%mud)
   u_drop = u_drop/max(1e0_CUSTOM_REAL,abs(T0(3,:))) 

   where (( f%theta>fi%idrop*f%dr ) .and. ( theta_old<fi%idrop*f%dr ) .and. (fi%idrop < NDROPS) )
     fi%mus    = fi%mud
     fi%mud    = f%mud - fi%idrop*(u_drop)      
     fi%theta  = fi%idrop*f%dr
     fi%idrop  = fi%idrop + 1
   end where
       
   mu    = fi%mus -((fi%mus-fi%mud)/f%dc) *(f%theta-fi%theta)
   mu    = max(mu,fi%mud)
!   mu    = max(mu,mumin) ! to avoid values lower than 0 .

  end function swf_multi_mu


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
  subroutine write_dataT_all(nt)

  integer, intent(in) :: nt
 
  integer :: i

  if (.not.allocated(faults)) return
  do i = 1,size(faults)
    call SCEC_write_dataT(faults(i)%dataT,faults(i)%dt,nt)
  enddo

  end subroutine write_dataT_all

!------------------------------------------------------------------------
  subroutine SCEC_write_dataT(dataT,DT,NT)

  type(dataT_type), intent(in) :: dataT
  real(kind=CUSTOM_REAL), intent(in) :: DT
  integer, intent(in) :: NT

  integer   :: i,k,IOUT
  character :: NTchar*5

  IOUT = 121 !WARNING: not very robust. Could instead look for an available ID

  write(NTchar,1) NT
  NTchar = adjustl(NTchar)

1 format(I5)  
 do i=1,dataT%npoin

      open(IOUT,file='OUTPUT_FILES/'//trim(dataT%name(i))//'.dat',status='replace')
      write(IOUT,*) "# problem=TPV15"
      write(IOUT,*) "# author=Galvez, Ampuero, Nissen-Meyer"
      write(IOUT,*) "# date=2011/xx/xx"
      write(IOUT,*) "# code=SPECFEM3D_FAULT "
      write(IOUT,*) "# code_version=1.1"
      write(IOUT,*) "# element_size=100 m  (*4 GLL nodes)"
      write(IOUT,*) "# time_step=",DT
      write(IOUT,*) "# num_time_steps=",NT
      write(IOUT,*) "# location=",trim(dataT%name(i))
      write(IOUT,*) "# Time series in 8 column of E15.7"
      write(IOUT,*) "# Column #1 = Time (s)"
      write(IOUT,*) "# Column #2 = horizontal right-lateral slip (m)"
      write(IOUT,*) "# Column #3 = horizontal right-lateral slip rate (m/s)"
      write(IOUT,*) "# Column #4 = horizontal right-lateral shear stress (MPa)"
      write(IOUT,*) "# Column #5 = vertical up-dip slip (m)"
      write(IOUT,*) "# Column #6 = vertical up-dip slip rate (m/s)"
      write(IOUT,*) "# Column #7 = vertical up-dip shear stress (MPa)"
      write(IOUT,*) "# Column #8 = normal stress (MPa)"
      write(IOUT,*) "#"
      write(IOUT,*) "# The line below lists the names of the data fields:"
      write(IOUT,*) "#t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress"
      write(IOUT,*) "#"
      do k=1,NT
        write(IOUT,'(8(E15.7))') k*DT, dataT%d1(k,i), dataT%v1(k,i), dataT%t1(k,i)/1.0e6_CUSTOM_REAL, &
                                         dataT%d2(k,i), dataT%v2(k,i), dataT%t2(k,i)/1.0e6_CUSTOM_REAL, &
                                         dataT%t3(k,i)/1.0e6_CUSTOM_REAL
      enddo
      close(IOUT)
  enddo

  end subroutine SCEC_write_dataT

!-------------------------------------------------------------------------------------------------

   subroutine SCEC_Write_RuptureTime(dataXZ,DT,NT,iflt)
 
  type(dataXZ_type), intent(in) :: dataXZ
  real(kind=CUSTOM_REAL), intent(in) :: DT
  integer, intent(in) :: NT,iflt
  
  integer   :: i,IOUT
  character(len=70) :: filename
    
  write(filename,"('OUTPUT_FILES/RuptureTime_Fault',I0)") iflt

  IOUT = 121 !WARNING: not very robust. Could instead look for an available ID
      
      open(IOUT,file=trim(filename),status='replace')
      write(IOUT,*) "# problem=TPV5"
      write(IOUT,*) "# author=Galvez, Ampuero, Tarje"
      write(IOUT,*) "# date=2011/xx/xx"
      write(IOUT,*) "# code=SPECFEM3D_FAULT"
      write(IOUT,*) "# code_version=1.1"
      write(IOUT,*) "# element_size=100 m  (*4 GLL nodes)"
      write(IOUT,*) "# time_step=",DT
      write(IOUT,*) "# num_time_steps=",NT
      write(IOUT,*) "# Column #1 = horizontal coordinate, distance along strike (m)"
      write(IOUT,*) "# Column #2 = vertical coordinate, distance down-dip (m)"
      write(IOUT,*) "# Column #3 = rupture time (s)"
      write(IOUT,*) "# x y z time"
     do i = 1,size(dataXZ%tRUP)
      write(IOUT,'(4(E15.7))') dataXZ%xcoord(i), dataXZ%ycoord(i), dataXZ%zcoord(i), dataXZ%tRUP(i)
     end do 

    close(IOUT)

   end subroutine SCEC_Write_RuptureTime

!-------------------------------------------------------------------------------------------------

  subroutine init_dataXZ(DataXZ,bc,nglob)

  type(dataXZ_type), intent(inout) :: DataXZ
  type(bc_dynflt_type) :: bc
  integer, intent(in) :: nglob

  allocate(DataXZ%stg(nglob))
  DataXZ%sta => bc%swf%theta
  DataXZ%d1 => bc%d(1,:)
  DataXZ%d2 => bc%d(2,:)
  DataXZ%v1 => bc%v(1,:)
  DataXZ%v2 => bc%v(2,:) 
  DataXZ%t1 => bc%t(1,:)
  DataXZ%t2 => bc%t(2,:)
  DataXZ%t3 => bc%t(3,:)
  DataXZ%xcoord => bc%coord(1,:) 
  DataXZ%ycoord => bc%coord(2,:)
  DataXZ%zcoord => bc%coord(3,:)
  allocate(DataXZ%tRUP(nglob))
  allocate(DataXZ%tPZ(nglob))

!Percy , setting up initial rupture time null for all faults.  
  DataXZ%tRUP = 0e0_CUSTOM_REAL
  DataXZ%tPZ  = 0e0_CUSTOM_REAL


  end subroutine init_dataXZ

!---------------------------------------------------------------
subroutine store_dataXZ(dataXZ,stg,dold,dnew,dc,vold,vnew,time,dt) 

  type(dataXZ_type), intent(inout) :: dataXZ
  real(kind=CUSTOM_REAL), dimension(:), intent(in) :: stg,dold,dnew,dc,vold,vnew
  real(kind=CUSTOM_REAL), intent(in) :: time,dt

  integer :: i

! "stg" : strength .
 
  dataXZ%stg   = stg

  do i = 1,size(stg)
   ! process zone time = first time when slip = dc  (break down process).
   ! with linear time interpolation
    if (dataXZ%tPZ(i)==0e0_CUSTOM_REAL) then
      if (dold(i)<=dc(i) .and. dnew(i) >= dc(i)) then
        dataXZ%tPZ(i) = time-dt*(dnew(i)-dc(i))/(dnew(i)-dold(i))
      endif
    endif
   ! rupture time = first time when slip velocity = vc
   ! with linear time interpolation
   ! vc should be pre-defined as input data .
  
    if (dataXZ%tRUP(i)==0e0_CUSTOM_REAL) then
      if (vold(i)<=V_RUPT .and. vnew(i)>=V_RUPT) dataXZ%tRUP(i)= time-dt*(vnew(i)-V_RUPT)/(vnew(i)-vold(i))
    endif
  enddo
  
! To do : add stress criteria (firs time strength is reached).

  ! note: the other arrays in dataXZ are pointers to arrays in bc
  !       they do not need to be updated here

  end subroutine store_dataXZ

!---------------------------------------------------------------
  subroutine write_dataXZ(dataXZ,itime,iflt)


  type(dataXZ_type), intent(in) :: dataXZ
  integer, intent(in) :: itime,iflt
   
  character(len=70) :: filename


  write(filename,"('OUTPUT_FILES/Snapshot',I0,'_F',I0,'.bin')") itime,iflt
! open(unit=IOUT, file= trim(filename), status='replace', form='formatted',action='write')
! NOTE : It had to be adopted formatted output to avoid conflicts readings with different 
!        compilers.

!  write(IOUT,"(5F24.15)") dataXZ%xcoord,dataXZ%ycoord,dataXZ%zcoord,dataXZ%v1,dataXZ%v2
 
  open(unit=IOUT, file= trim(filename), status='replace', form='unformatted',action='write')

  write(IOUT) dataXZ%xcoord
  write(IOUT) dataXZ%ycoord
  write(IOUT) dataXZ%zcoord
  write(IOUT) dataXZ%d1
  write(IOUT) dataXZ%d2
  write(IOUT) dataXZ%v1
  write(IOUT) dataXZ%v2
  write(IOUT) dataXZ%t1
  write(IOUT) dataXZ%t2
  write(IOUT) dataXZ%t3
  write(IOUT) dataXZ%sta
  write(IOUT) dataXZ%stg
  write(IOUT) dataXZ%tRUP
  write(IOUT) dataXZ%tPZ
  close(IOUT)

  end subroutine write_dataXZ

!---------------------------------------------------------------
  subroutine write_init_fault_par(bc,iflt)


  type(bc_dynflt_type), intent(inout) :: bc
  integer, intent(in) :: iflt
   
  character(len=70) :: filename


  write(filename,"('OUTPUT_FILES/Fault_init_par','_F',I0,'.bin')") iflt
! open(unit=IOUT, file= trim(filename), status='replace', form='formatted',action='write')
! NOTE : It had to be adopted formatted output to avoid conflicts readings with different 
!        compilers.

!  write(IOUT,"(5F24.15)") dataXZ%xcoord,dataXZ%ycoord,dataXZ%zcoord,dataXZ%v1,dataXZ%v2
 
  open(unit=IOUT, file= trim(filename), status='replace', form='unformatted',action='write')

  write(IOUT) bc%coord(1,:)
  write(IOUT) bc%coord(2,:)
  write(IOUT) bc%coord(3,:)
  write(IOUT) bc%swf%mus 
  write(IOUT) bc%swf%mud
  write(IOUT) bc%swf%Dc
  write(IOUT) bc%swf%Dr
  write(IOUT) bc%T0(1,:)
  write(IOUT) bc%T0(2,:)
  write(IOUT) bc%T0(3,:)

  close(IOUT)

  end subroutine write_init_fault_par 

!---

end module fault_solver
