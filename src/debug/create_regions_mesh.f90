!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!
!-------------------------------------------------------------------------------------------------
!

! main routine

subroutine create_regions_mesh_ext(ibool, &
                        xstore,ystore,zstore,nspec,npointot,myrank,LOCAL_PATH, &
                        nnodes_ext_mesh,nelmnts_ext_mesh, &
                        nodes_coords_ext_mesh, elmnts_ext_mesh, &
                        max_static_memory_size, mat_ext_mesh, materials_ext_mesh, &
                        nmat_ext_mesh, undef_mat_prop, nundefMat_ext_mesh, &
                        num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh, my_nelmnts_neighbours_ext_mesh, &
                        my_interfaces_ext_mesh, &
                        ibool_interfaces_ext_mesh, nibool_interfaces_ext_mesh, &
                        nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, &
                        NSPEC2D_BOTTOM, NSPEC2D_TOP,&
                        ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
                        nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax,&
                        nodes_ibelm_bottom,nodes_ibelm_top, &
                        SAVE_MESH_FILES,nglob, &
                        ANISOTROPY,NPROC,OCEANS, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)

! create the different regions of the mesh

  use create_regions_mesh_ext_par
  use fault_object, only: fault_read_input,fault_setup, &
                          fault_save_arrays,fault_save_arrays_test,&
                          nnodes_coords_open,nodes_coords_open,ANY_FAULT_IN_THIS_PROC,&
                          ANY_FAULT 

  implicit none

! number of spectral elements in each block
  integer :: nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer :: npointot

! proc numbers for MPI
  integer :: myrank
  integer :: NPROC

  character(len=256) :: LOCAL_PATH

!  data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! static memory size needed by the solver
  double precision :: max_static_memory_size

  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh

! material properties
  integer :: nmat_ext_mesh,nundefMat_ext_mesh
  double precision, dimension(6,nmat_ext_mesh) :: materials_ext_mesh
  character (len=30), dimension(6,nundefMat_ext_mesh):: undef_mat_prop

!  double precision, external :: materials_ext_mesh

! MPI communication
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: my_interfaces_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh

! absorbing boundaries
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
  integer, dimension(nspec2D_xmin)  :: ibelm_xmin
  integer, dimension(nspec2D_xmax)  :: ibelm_xmax
  integer, dimension(nspec2D_ymin)  :: ibelm_ymin
  integer, dimension(nspec2D_ymax)  :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM)  :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP)  :: ibelm_top
  ! node indices of boundary faces
  integer, dimension(4,nspec2D_xmin)  :: nodes_ibelm_xmin
  integer, dimension(4,nspec2D_xmax)  :: nodes_ibelm_xmax
  integer, dimension(4,nspec2D_ymin)  :: nodes_ibelm_ymin
  integer, dimension(4,nspec2D_ymax)  :: nodes_ibelm_ymax
  integer, dimension(4,NSPEC2D_BOTTOM)  :: nodes_ibelm_bottom
  integer, dimension(4,NSPEC2D_TOP)  :: nodes_ibelm_top

  integer :: nglob

  logical :: SAVE_MESH_FILES
  logical :: ANISOTROPY
  logical :: OCEANS

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

! local parameters
! static memory size needed by the solver
  double precision :: static_memory_size
  real(kind=CUSTOM_REAL) :: model_speed_max

! for vtk output
!  character(len=256) prname_file
!  integer,dimension(:),allocatable :: itest_flag
!  integer, dimension(:), allocatable :: elem_flag

! For Piero Basini :
! integer :: doubling_value_found_for_Piero
!   double precision :: xmesh,ymesh,zmesh
!   double precision :: rho,vp,vs

!   integer,dimension(nspec) ::  idoubling
!   integer :: doubling_value_found_for_Piero
!   integer, parameter :: NUMBER_OF_STATIONS = 6
!   double precision, parameter :: RADIUS_TO_EXCLUDE = 250.d0
!   double precision, dimension(NUMBER_OF_STATIONS) :: utm_x_station,utm_y_station

!   logical :: is_around_a_station
!   integer :: istation

! ! store bedrock values
!   integer ::  icornerlat,icornerlong
!   double precision ::  lat,long,elevation_bedrock
!   double precision ::  lat_corner,long_corner,ratio_xi,ratio_eta
!real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ibedrock


! initializes arrays
  call sync_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...allocating arrays '
  endif
  call crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                        nspec2D_bottom,nspec2D_top,ANISOTROPY)

 ! if faults exist this reads nodes_coords_open
  call fault_read_input(prname,NDIM)

  call sync_all()
  if (myrank == 0) write(IMAIN,*) '  ...setting up jacobian '
  if (ANY_FAULT_IN_THIS_PROC) then
   ! compute jacobians with fault open and *store needed for ibool.
    call crm_ext_setup_jacobian(myrank, &
                         xstore,ystore,zstore,nspec, &
                         nodes_coords_open, nnodes_coords_open,&
                         elmnts_ext_mesh,nelmnts_ext_mesh)
  else ! with no fault
    call crm_ext_setup_jacobian(myrank, &
                        xstore,ystore,zstore,nspec, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,&
                        elmnts_ext_mesh,nelmnts_ext_mesh)
  endif

 ! create ibool with faults open
  call sync_all()
  if (myrank == 0) write(IMAIN,*) '  ...indexing global points'
  if (ANY_FAULT_IN_THIS_PROC) then
    call crm_ext_setup_indexing(ibool, &
                       xstore,ystore,zstore,nspec,nglob,npointot, &
                       nnodes_coords_open,nodes_coords_open,myrank)
  else ! with no fault
    call crm_ext_setup_indexing(ibool, &
                      xstore,ystore,zstore,nspec,nglob,npointot, &
                      nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)
  end if

  if (ANY_FAULT) then
   ! recalculate *store with faults closed 
    call sync_all()
    if (myrank == 0) write(IMAIN,*) '  ... resetting up jacobian in fault domains'
    if (ANY_FAULT_IN_THIS_PROC) call crm_ext_setup_jacobian(myrank, &
                           xstore,ystore,zstore,nspec, &
                           nodes_coords_ext_mesh,nnodes_ext_mesh,&
                           elmnts_ext_mesh,nelmnts_ext_mesh)
   ! at this point (xyz)store_dummy are still open
    call fault_setup (ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                    xstore,ystore,zstore,nspec,nglob,myrank)
   ! this closes (xyz)store_dummy
  endif


! sets up MPI interfaces between partitions
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...preparing MPI interfaces '
  endif
  call get_MPI(myrank,nglob,nspec,ibool, &
                        nelmnts_ext_mesh,elmnts_ext_mesh, &
                        my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                        ibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh,&
                        my_neighbours_ext_mesh,NPROC)

! sets material velocities
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...determining velocity model'
  endif
  call get_model(myrank,nspec,ibool,mat_ext_mesh,nelmnts_ext_mesh, &
                        materials_ext_mesh,nmat_ext_mesh, &
                        undef_mat_prop,nundefMat_ext_mesh, &
                        ANISOTROPY)

! sets up absorbing/free surface boundaries
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up absorbing boundaries '
  endif
  call get_absorbing_boundary(myrank,nspec,nglob,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            nspec2D_bottom,nspec2D_top)

! sets up acoustic-elastic coupling surfaces
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...detecting acoustic-elastic surfaces '
  endif
  call get_coupling_surfaces(myrank, &
                        nspec,nglob,ibool,NPROC, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh)


! creates mass matrix
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating mass matrix '
  endif
  call create_mass_matrices(nglob,nspec,ibool)

! creates ocean load mass matrix
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating ocean load mass matrix '
  endif
  call create_mass_matrices_ocean_load(nglob,nspec,ibool,OCEANS,&
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)

! saves the binary mesh files
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...saving databases'
  endif
  !call create_name_database(prname,myrank,LOCAL_PATH)
  call save_arrays_solver_ext_mesh(nspec,nglob, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore, &
                        jacobianstore, rho_vp,rho_vs,iflag_attenuation_store, &
                        rhostore,kappastore,mustore, &
                        rmass,rmass_acoustic,rmass_solid_poroelastic,rmass_fluid_poroelastic, &
                        OCEANS,rmass_ocean_load,NGLOB_OCEAN,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec,num_abs_boundary_faces, &
                        free_surface_normal,free_surface_jacobian2Dw, &
                        free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
                        coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                        coupling_ac_el_ijk,coupling_ac_el_ispec,num_coupling_ac_el_faces, &
                        num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                        max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                        prname,SAVE_MESH_FILES,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  call fault_save_arrays_test(prname,IOUT)  ! for debugging
  call fault_save_arrays(prname,IOUT)

! computes the approximate amount of static memory needed to run the solver
  call memory_eval(nspec,nglob,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh,static_memory_size)
  call max_all_dp(static_memory_size, max_static_memory_size)

! checks the mesh, stability and resolved period
  call sync_all()
  call check_mesh_resolution(myrank,nspec,nglob,ibool,&
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            kappastore,mustore,rho_vp,rho_vs, &
                            -1.0d0, model_speed_max )

! VTK file output
!  if( SAVE_MESH_FILES ) then
!    ! saves material flag assigned for each spectral element into a vtk file
!    prname_file = prname(1:len_trim(prname))//'material_flag'
!    allocate(elem_flag(nspec))
!    elem_flag(:) = mat_ext_mesh(1,:)
!    call write_VTK_data_elem_i(nspec,nglob, &
!            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!            elem_flag,prname_file)
!    deallocate(elem_flag)
!
!    !plotting abs boundaries
!    !  allocate(itest_flag(nspec))
!    !  itest_flag(:) = 0
!    !  do ispec=1,nspec
!    !    if( iboun(1,ispec) ) itest_flag(ispec) = 1
!    !  enddo
!    !  prname_file = prname(1:len_trim(prname))//'iboundary1_flag'
!    !  call write_VTK_data_elem_i(nspec,nglob, &
!    !            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!    !            itest_flag,prname_file)
!    !  deallocate(itest_flag)
!  endif

! AVS/DX file output
! create AVS or DX mesh data for the slice, edges and faces
!  if(SAVE_MESH_FILES) then
! check: no idoubling
!    call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
!    call write_AVS_DX_mesh_quality_data(prname,nspec,xstore,ystore,zstore, &
!                   kappastore,mustore,rhostore)
! check: no iMPIcut_xi,iMPIcut_eta,idoubling
!    call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
!              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
! check: no idoubling
!    call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
!              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
!  endif

! cleanup
  if( .not. SAVE_MOHO_MESH ) deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  deallocate(xixstore,xiystore,xizstore,&
              etaxstore,etaystore,etazstore,&
              gammaxstore,gammaystore,gammazstore)
  deallocate(jacobianstore,iflag_attenuation_store)
  deallocate(kappastore,mustore,rho_vp,rho_vs)

end subroutine create_regions_mesh_ext

!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                        nspec2D_bottom,nspec2D_top,ANISOTROPY)

  use create_regions_mesh_ext_par
  implicit none

  integer :: nspec,myrank
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            nspec2D_bottom,nspec2D_top

  character(len=256) :: LOCAL_PATH

  logical :: ANISOTROPY

! local parameters
  integer :: ier

! memory test
!  logical,dimension(:),allocatable :: test_mem
!
! tests memory availability (including some small buffer of 10*1024 byte)
!  allocate( test_mem(int(max_static_memory_size)+10*1024),stat=ier)
!  if(ier /= 0) then
!    write(IMAIN,*) 'error: try to increase the available process stack size by'
!    write(IMAIN,*) '       ulimit -s **** '
!    call exit_MPI(myrank,'not enough memory to allocate arrays')
!  endif
!  test_mem(:) = .true.
!  deallocate( test_mem, stat=ier)
!  if(ier /= 0) call exit_MPI(myrank,'error to allocate arrays')
!  call sync_all()

  allocate( xelm(NGNOD),yelm(NGNOD),zelm(NGNOD),stat=ier)

  allocate( iflag_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ))

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ))

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), &
          dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)

! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ), &
          shape2D_y(NGNOD2D,NGLLX,NGLLZ), &
          shape2D_bottom(NGNOD2D,NGLLX,NGLLY), &
          shape2D_top(NGNOD2D,NGLLX,NGLLY), stat=ier)

  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ), &
          dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
          dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY), &
          dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)

  allocate(wgllwgll_xy(NGLLX,NGLLY), &
          wgllwgll_xz(NGLLX,NGLLZ), &
          wgllwgll_yz(NGLLY,NGLLZ),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec), &
          rho_vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
          kappastore(NGLLX,NGLLY,NGLLZ,nspec), &
          mustore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
          !vpstore(NGLLX,NGLLY,NGLLZ,nspec), &
          !vsstore(NGLLX,NGLLY,NGLLZ,nspec),
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! arrays with mesh parameters
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec), &
          xiystore(NGLLX,NGLLY,NGLLZ,nspec), &
          xizstore(NGLLX,NGLLY,NGLLZ,nspec), &
          etaxstore(NGLLX,NGLLY,NGLLZ,nspec), &
          etaystore(NGLLX,NGLLY,NGLLZ,nspec), &
          etazstore(NGLLX,NGLLY,NGLLZ,nspec), &
          gammaxstore(NGLLX,NGLLY,NGLLZ,nspec), &
          gammaystore(NGLLX,NGLLY,NGLLZ,nspec), &
          gammazstore(NGLLX,NGLLY,NGLLZ,nspec), &
          jacobianstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! absorbing boundary
  ! absorbing faces
  num_abs_boundary_faces = nspec2D_xmin + nspec2D_xmax + nspec2D_ymin + nspec2D_ymax + nspec2D_bottom
  ! adds faces of free surface if it also absorbs
  if( ABSORB_FREE_SURFACE ) num_abs_boundary_faces = num_abs_boundary_faces + nspec2D_top

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate( abs_boundary_ispec(num_abs_boundary_faces), &
           abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! free surface
  ! free surface faces
  if( ABSORB_FREE_SURFACE ) then
    ! no free surface - uses a dummy size
    num_free_surface_faces = 1
  else
    num_free_surface_faces = nspec2D_top
  endif

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate( free_surface_ispec(num_free_surface_faces), &
           free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces), &
           free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces), &
           free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with anisotropy
  if( ANISOTROPY ) then
    NSPEC_ANISO = nspec
  else
    NSPEC_ANISO = 1
  endif
  allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! material flags
  allocate( ispec_is_acoustic(nspec), &
           ispec_is_elastic(nspec), &
           ispec_is_poroelastic(nspec), stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

end subroutine crm_ext_allocate_arrays


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_setup_jacobian(myrank, &
                        xstore,ystore,zstore,nspec, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,&
                        elmnts_ext_mesh,nelmnts_ext_mesh)

  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: nspec

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! proc numbers for MPI
  integer :: myrank

! local parameters
  integer :: ispec,ia,i,j,k

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY)

! 2D weights
  do j=1,NGLLY
    do i=1,NGLLX
      wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo
  do k=1,NGLLZ
    do i=1,NGLLX
      wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo
  do k=1,NGLLZ
    do j=1,NGLLY
      wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

! point locations
  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0
  do ispec = 1, nspec
    !call get_xyzelm(xelm, yelm, zelm, ispec, elmnts_ext_mesh, nodes_coords_ext_mesh, nspec, nnodes_ext_mesh)
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
    enddo
    ! CUBIT should provide a mesh ordering such that the 3D jacobian is defined
    ! (otherwise mesh would be degenerated)
    call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
                      etaxstore,etaystore,etazstore, &
                      gammaxstore,gammaystore,gammazstore,jacobianstore, &
                      xstore,ystore,zstore, &
                      xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)
  enddo

end subroutine crm_ext_setup_jacobian


!
!-------------------------------------------------------------------------------------------------

subroutine crm_ext_setup_indexing(ibool, &
                            xstore,ystore,zstore,nspec,nglob,npointot, &
                            nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)

! creates global indexing array ibool

  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: nspec,nglob,npointot,myrank

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

! data from the external mesh
  integer :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh

! local parameters
! variables for creating array ibool
  double precision, dimension(:), allocatable :: xp,yp,zp
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg

  integer :: ieoff,ilocnum,ier
  integer :: i,j,k,ispec,iglobnum

! allocate memory for arrays
  allocate(locval(npointot), &
          ifseg(npointot), &
          xp(npointot), &
          yp(npointot), &
          zp(npointot),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! creates temporary global point arrays
  locval = 0
  ifseg = .false.
  xp = 0.d0
  yp = 0.d0
  zp = 0.d0

  do ispec=1,nspec
    ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
    ilocnum = 0
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ilocnum = ilocnum + 1
          xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
          yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
          zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

! gets ibool indexing from local (GLL points) to global points
  call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot, &
       minval(nodes_coords_ext_mesh(1,:)),maxval(nodes_coords_ext_mesh(1,:)))

!- we can create a new indirect addressing to reduce cache misses
  call get_global_indirect_addressing(nspec,nglob,ibool)

!cleanup
  deallocate(xp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(locval,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if(ier /= 0) stop 'error in deallocate'

! unique global point locations
  allocate(xstore_dummy(nglob), &
          ystore_dummy(nglob), &
          zstore_dummy(nglob),stat=ier)
  if(ier /= 0) stop 'error in allocate'
  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              xstore_dummy(iglobnum) = xstore(i,j,k,ispec)
              ystore_dummy(iglobnum) = ystore(i,j,k,ispec)
              zstore_dummy(iglobnum) = zstore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

end subroutine crm_ext_setup_indexing
