module fault_scotch

  implicit none

  private 
  type fault_type
    private 
    integer :: nspec
    integer, dimension(:), pointer  :: ispec1, ispec2 !, iface1, iface2 
    integer, dimension(:,:), pointer  :: inodes1, inodes2 
  end type fault_type

  type(fault_type), allocatable, save :: faults(:) 
  double precision, dimension(:,:), allocatable, save :: nodes_coords_open
 
 
  integer, parameter :: long = SELECTED_INT_KIND(18)

  double precision, parameter :: FAULT_GAP_TOLERANCE = 1.0d0

  public :: read_fault_files, fault_collect_elements, close_faults, write_fault_database, &
            save_nodes_coords, nodes_coords_open, faults

CONTAINS 
!==========================================================================================

  Subroutine read_fault_files(localpath_name)

  character(len=256),intent(in) :: localpath_name    
  integer :: nbfaults, iflt, ier 

  open(101,file='../DATA/FAULT/Par_file_faults.in',status='old',action='read',iostat=ier)
  if (ier==0) then 
    read(101,*) nbfaults
  else
    nbfaults = 0
    print *, 'Par_file.in not found: assume no faults'
  endif
  close(101)

  if (nbfaults>0) then
    allocate(faults(nbfaults))
    do iflt = 1 , nbfaults 
     call read_single_fault_file(faults(iflt),iflt,localpath_name)
    enddo
  endif

  end subroutine read_fault_files


!---------------------------------------------------------------------------------------------------

  Subroutine read_single_fault_file(f,ifault,localpath_name)

  type(fault_type), intent(inout) :: f
  character(len=256),intent(in) :: localpath_name 
 
  character(len=256) :: filename  
  integer,intent(in) :: ifault
  character(len=5) :: NTchar
  integer :: e,ier,nspec_side1,nspec_side2 
  
  write(NTchar,'(I5)') ifault
  NTchar = adjustl(NTchar)
   
  filename = localpath_name(1:len_trim(localpath_name))//'/fault_file_'//&
             NTchar(1:len_trim(NTchar))//'.dat'
  filename = adjustl(filename)
  ! reads fault elements and nodes
  open(unit=101, file=filename, status='old', form='formatted', iostat = ier)
  if( ier /= 0 ) then
    write(6,*) 'Fatal error: file '//filename//' not found' 
    write(6,*) 'Abort'
    stop
  endif
       
  read(101,*) nspec_side1,nspec_side2
  if (nspec_side1 /= nspec_side2) stop 'Number of fault nodes at do not match.'
  f%nspec = nspec_side1
  allocate(f%ispec1(f%nspec))
  allocate(f%ispec2(f%nspec))
  allocate(f%inodes1(4,f%nspec))
  allocate(f%inodes2(4,f%nspec))
 ! format: 
 ! number of elements side at side 1 and side 2 
 ! #id_(element containing the face) #id_node1_face .. #id_node4_face
 ! First for all faces on side 1, then side 2
 ! Note: element ids start at 1, not 0. Check that in cubit2specfem3d.py
  do e=1,f%nspec
    read(101,*) f%ispec1(e),f%inodes1(:,e)
!      print*, f%ispec1(e),f%inodes1(:,e) !TEST
  enddo
  do e=1,f%nspec
    read(101,*) f%ispec2(e),f%inodes2(:,e)
!      print*, f%ispec2(e),f%inodes2(:,e) !TEST
  enddo
 ! If we ever figure out how to export "ifaces" from CUBIT:
  !allocate(f%iface1(f%nspec))
  !allocate(f%iface2(f%nspec))
  !do e=1,f%nspec
  !  read(101,*) f%ispec1(e),f%ispec2(e),f%iface1(e),f%iface2(e)
  !enddo

  close(101)

  end Subroutine read_single_fault_file


! ---------------------------------------------------------------------------------------------------
! Saving nodes_coords to be used in SESAME for ibool_fault_side1 and side2
   subroutine save_nodes_coords(nodes_coords,nnodes)
   
   integer, intent(in) :: nnodes
   double precision, dimension(3,nnodes), intent(in) :: nodes_coords
  
   allocate(nodes_coords_open(3,nnodes)) 
   nodes_coords_open = nodes_coords 
 
   end subroutine save_nodes_coords   

! ---------------------------------------------------------------------------------------------------

  subroutine close_faults(nodes_coords,elmnts,nelmnts,nnodes,esize)
    
  integer ,intent(in) :: nnodes, esize
  integer(long), intent(in) :: nelmnts
  double precision, dimension(3,nnodes), intent(inout) :: nodes_coords
  integer, dimension(esize,nelmnts), intent(in) :: elmnts

  integer  :: iflt

  do iflt=1,size(faults)
    call close_fault_single(faults(iflt)%ispec1,faults(iflt)%ispec2,&
                            elmnts,nodes_coords,nnodes,esize,nelmnts)
  enddo

  end subroutine close_faults

! ---------------------------------------------------------------------------------------------------
!jpa: to do this much faster:
!     1. create a list of unique nodes from inodes1 and inodes2, 
!          inodes1_u = unique(isort1)
!          inodes2_u = unique(isort2)
!     2. sort the nodes by coordinates. Now both faces correspond.
!          [coord,k1] = sort(coords(inodes1_u))
!          k1 = inodes1_u(k1)
!          [coord,k2] = sort(coords(inodes2_u))
!          k2 = inodes2_u(k2)
!     3. set the coordinates on both sides equal to their average
!          coords(k1) = 0.5*( coords(k1)+coords(k2) )
!          coords(k2) = coords(k1)
  subroutine close_fault_single(ispec1,ispec2,elmnts,nodes_coords,nnodes,esize,nelmnts)
 
  integer ,intent(in)  :: nnodes, esize
  integer(long), intent(in) :: nelmnts
  integer, dimension(esize,nelmnts), intent(in) :: elmnts
  integer , dimension(:), intent(in) :: ispec1,ispec2
  double precision,dimension(3,nnodes), intent(inout) :: nodes_coords 
    
  double precision, dimension(3) :: xyz_1, xyz_2 
  double precision, dimension(3) :: xyz
  
  double precision :: dist
  integer :: iglob1, iglob2, i, j, k1, k2
  logical :: found_it

  do i = 1,size(ispec2)
    do k2=1,esize
      iglob2 = elmnts(k2,ispec2(i))
      found_it = .false.
      xyz_2 = nodes_coords(:,iglob2)

      do j = 1,size(ispec1)
        do k1=1,esize
          iglob1 = elmnts(k1,ispec1(j))
          xyz_1 = nodes_coords(:,iglob1)
          xyz = xyz_2-xyz_1
          dist = xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)
          dist = sqrt(dist)
  
       !jpa: Closing nodes that are already closed is not a problem
       !jpa: I process them again to leave the loop as early as possible
       !jpa: and to test if a facing node was found (See below). 
   
          if (dist <= FAULT_GAP_TOLERANCE) then
            xyz =  (xyz_1 + xyz_2)*0.5d0
            nodes_coords(:,iglob2) = xyz
            nodes_coords(:,iglob1) = xyz
            found_it = .true.
            exit 
          endif
 
        enddo
        if (found_it) exit 
      enddo

    enddo
       ! if (.not.found_it) then
       ! If the fault is very complicated (non-planar) the meshes of the two fault sides 
       ! can be very unstructured and might not match (unless you enforce it in CUBIT). 
       ! That is unlikely because the fault gap is tiny, but we never know.
       !   stop 'Inconsistent fault mesh: corresponding node in the other fault face was not found'
       ! endif
  enddo
  
  end subroutine close_fault_single

! ---------------------------------------------------------------------------------------------------
  !--------------------------------------------------
  ! Repartitioning : two coupled faultside1/side2 elements are transfered to the same partition
  !--------------------------------------------------

 Subroutine fault_collect_elements(nelmnts,nnodes,elmnts, &
                                   sup_neighbour,esize,nsize,nproc,part)

! INPUTS
  integer(long), intent(in) :: nelmnts,nsize,sup_neighbour 
  integer, dimension(0:esize*nelmnts-1),intent(in)  :: elmnts
  integer, intent(in)  :: nnodes, nproc, esize
! OUTPUTS :
  integer, dimension(0:nelmnts-1),intent(inout)    :: part
! VARIABLES:
  logical, dimension(nelmnts)  :: is_on_fault
  integer :: iflt
  
  is_on_fault = .false.
  do iflt=1,size(faults)
    is_on_fault(faults(iflt)%ispec1) = .true.
    is_on_fault(faults(iflt)%ispec2) = .true.
  end do
  call fault_repartition (nelmnts, nnodes, elmnts, sup_neighbour, nsize, &
                          nproc, part, is_on_fault,esize)

  end Subroutine fault_collect_elements

! ---------------------------------------------------------------------------------------------------

  Subroutine fault_repartition (nelmnts, nnodes, elmnts, sup_neighbour, nsize, &
                        nproc, part, is_on_fault, esize)

!  INDIVIDUAL FAULT REPARTITION

!     part : iproc number of processor partioned. It will altered patching fault elements into the same partion.  
!     Part, once is altered , will be input for write_partition_database.

!INPUTS
  integer(long), intent(in) :: nelmnts,sup_neighbour,nsize
  integer, intent(in)  :: nnodes, nproc, esize 
  integer, dimension(0:esize*nelmnts-1), intent(in) :: elmnts
  logical , dimension(nelmnts), intent(in) :: is_on_fault
!OUTPUTS :
  integer, dimension(0:nelmnts-1), intent(inout)    :: part

!LOCAL VARIABLES :
  integer, dimension(0:nelmnts)                  :: xadj
  integer, dimension(0:sup_neighbour*nelmnts-1)  :: adjncy
  integer, dimension(0:nnodes-1)                 :: nnodes_elmnts
  integer, dimension(0:nsize*nnodes-1)           :: nodes_elmnts
  integer  :: max_neighbour       

!SHILDING 
  integer  :: i,j, ipart,nproc_null,nproc_null_final
  integer  :: el, el_1, el_2, k1, k2, k,e,iflt,inode
  logical  :: is_repartitioned
  integer, dimension(:), allocatable :: elem_proc_null

 ! downloading processor 0
  nproc_null = count( part == 0 )

  print*, 'Elements proc = 0 redistributed in [{nproc}- nproc0] :'
  print*, nproc_null

  if ( nproc_null /= 0 ) then
 
    allocate(elem_proc_null(nproc_null))
   ! Filling up proc = 0 elements
    nproc_null = 0
    do i = 0,nelmnts-1
      if ( part(i) == 0 ) then
        nproc_null = nproc_null + 1
        elem_proc_null(nproc_null) = i
      end if
    end do     
   ! Redistributing proc-0 elements on the rest of processors
   !jpa: why do this? does it always help balancing ?
   !pgb: Yes, bulk elements in processor 0 are taken out and redistributed.
   !pgb: leaving more space for fault elements. 
   !jpa: But if the number of fault elements is much smaller than nproc_null
   !     we will end up with a very UNbalanced proc 0 !
   !pgb: Solution , fault parellelization . 
    ipart=1
    if (nproc > 1) then 
      do i = 1, nproc_null
        part(elem_proc_null(i)) = ipart
        if ( ipart == nproc-1 ) ipart = 0
        ipart = ipart +1
      end do
    end if
    deallocate(elem_proc_null)
  endif
!  call mesh2dual_ncommonnodes_fault(nelmnts, nnodes, nsize, sup_neighbour, &
!                                elmnts, xadj, adjncy, nnodes_elmnts, &
!                                nodes_elmnts, max_neighbour, 4, esize)
  
    ! coupled elements
    !  ---------------
    ! Allocating neighbours with shared fault faces.
    ! This is done in order to not break MPI interfaces.
    ! Dirty implementation.
    ! Fault faces are shield by double coarse neighbours.
    !              
    !                 1   2
    !             1   1   2   2
    !          1  1  [1] [2]  2  2
    !             1  1    2   2      
    !                1    2   
                 
    ! Allocating elements with double shield layer
!
!  ===========   FAULT SHIELD double-layer =========== !
!  print *, "Fault shield double-layer :"
!  do el = 0, nelmnts-1
!    if ( is_on_fault(el+1) ) then
!      part(el) = 0
!      do k1 = xadj(el), xadj(el+1) - 1
!        el_1 = adjncy(k1) 
!        part(el_1) = 0
!        do k2 = xadj(el_1), xadj(el_1+1) - 1
!          el_2 = adjncy(k2) 
!          part(el_2) = 0
!        enddo
!      enddo
!    endif
!  enddo
! ===================================================== !

! Fault zone layer = the set of elements that contain at least one fault node
  print *, "Fault zone layer :"

! List of elements per node
!  nnodes_elmnts(i) = number of elements containing node #i (i=0:nnodes-1)     
!  nodes_elmnts(nsize*i:nsize*i+nnodes_elmnts(i)-1) = index of elements (starting at 0) containing node #i  
!  nsize = maximun number of elements in a node.
!  esize = nodes per element. 

  nnodes_elmnts(:) = 0
  nodes_elmnts(:)  = 0

  do i = 0, esize*nelmnts-1
    nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
    nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
  end do

  do iflt=1,size(faults)
    do e=1,faults(iflt)%nspec
      do k=1,4

        inode = faults(iflt)%inodes1(k,e)-1  ! node index, starting at 0
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        part( nodes_elmnts(k1:k2) ) = 0
        inode = faults(iflt)%inodes2(k,e)-1 
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        part( nodes_elmnts(k1:k2) ) = 0

      end do
    end do
  end do

  nproc_null_final = count( part == 0 )
  print *, nproc_null_final 

  end subroutine fault_repartition

! ---------------------------------------------------------------------------------------------------
! See subroutine write_boundaries_database in part_decompose_mesh_SCOTCH.f90
!
! File format: 
! one block for each fault
! first line of each block = number of fault elements in this processor
! next lines: #id_(element containing the face) #id_node1_face .. #id_node4_face
! first for all faces on side 1, then side 2

  subroutine write_fault_database(IIN_database, iproc, nelmnts, glob2loc_elmnts, &
                          glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, part)

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc
  integer(long), intent(in) :: nelmnts
  integer, dimension(0:nelmnts-1), intent(in)  :: part
  integer, dimension(0:nelmnts-1), intent(in)  :: glob2loc_elmnts
  integer, dimension(0:), intent(in) :: glob2loc_nodes_nparts
  integer, dimension(0:), intent(in) :: glob2loc_nodes_parts
  integer, dimension(0:), intent(in) :: glob2loc_nodes

  integer  :: i,j,k,iflt,e
  integer  :: nspec_fault_1,nspec_fault_2
  integer :: loc_nodes(4),inodes(4)

  do iflt=1,size(faults)
 
   ! get number of fault elements in this partition
    nspec_fault_1 = count( part(faults(iflt)%ispec1-1) == iproc )
    nspec_fault_2 = count( part(faults(iflt)%ispec2-1) == iproc )
    if (nspec_fault_1 /= nspec_fault_2) then
      print *, 'Fault # ',iflt,', proc # ',iproc
      print *, '  ispec1 : ', nspec_fault_1
      print *, '  ispec2 : ', nspec_fault_2
      print *, 'Fatal error: Number of fault elements do not coincide. Abort.'
      stop 
    end if
    write(IIN_database,*) nspec_fault_1

   ! if no fault element in this partition, move to next fault
    if (nspec_fault_1==0) cycle 

   ! export fault element data, side 1
    do i=1,faults(iflt)%nspec
      e = faults(iflt)%ispec1(i)
      if (part(e-1) == iproc) then
        inodes = faults(iflt)%inodes1(:,i)
        do k=1,4
          do j = glob2loc_nodes_nparts(inodes(k)-1), glob2loc_nodes_nparts(inodes(k))-1
            if (glob2loc_nodes_parts(j) == iproc ) then
              loc_nodes(k) = glob2loc_nodes(j) + 1
            end if
          end do
        end do
        write(IIN_database,*) glob2loc_elmnts(e-1)+1, loc_nodes
      end if
    enddo

   ! export fault element data, side 2
    do i=1,faults(iflt)%nspec
      e = faults(iflt)%ispec2(i)
      if(part(e-1) == iproc) then
        inodes = faults(iflt)%inodes2(:,i)
        do k=1,4
          do j = glob2loc_nodes_nparts(inodes(k)-1), glob2loc_nodes_nparts(inodes(k))-1
            if (glob2loc_nodes_parts(j) == iproc ) then
              loc_nodes(k) = glob2loc_nodes(j)+1
            end if
          end do
        end do
        write(IIN_database,*) glob2loc_elmnts(e-1)+1, loc_nodes
      end if
    enddo

  enddo

  end subroutine write_fault_database


! ---------------------------------------------------------------------------------------------------
! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
! Taken from part_decomposition_SCOTCH.f90 routine.
!----------------------------------------------------------------------------------------------------
  subroutine mesh2dual_ncommonnodes_fault(nelmnts, nnodes, nsize, sup_neighbour, elmnts,&
                        xadj, adjncy, &
                        nnodes_elmnts, nodes_elmnts, &
                        max_neighbour, ncommonnodes,esize)

    integer(long), intent(in)  :: nelmnts
    integer, intent(in)  :: nnodes
    integer(long), intent(in)  :: nsize
    integer(long), intent(in)  :: sup_neighbour
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts    
    
    integer, dimension(0:nelmnts)  :: xadj
    integer, dimension(0:sup_neighbour*nelmnts-1)  :: adjncy
    integer, dimension(0:nnodes-1)  :: nnodes_elmnts
    integer, dimension(0:nsize*nnodes-1)  :: nodes_elmnts
    integer, intent(out) :: max_neighbour
    integer, intent(in)  :: ncommonnodes,esize

    ! local parameters
    integer  :: i, j, k, l, m, nb_edges
    logical  ::  is_neighbour
    integer  :: num_node, n
    integer  :: elem_base, elem_target
    integer  :: connectivity


    ! initializes
    xadj(:) = 0
    adjncy(:) = 0
    nnodes_elmnts(:) = 0
    nodes_elmnts(:) = 0
    nb_edges = 0

    ! list of elements per node
    do i = 0, esize*nelmnts-1
       nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
       nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
    end do

    ! checking which elements are neighbours ('ncommonnodes' criteria)
    do j = 0, nnodes-1
       do k = 0, nnodes_elmnts(j)-1
          do l = k+1, nnodes_elmnts(j)-1

             connectivity = 0
             elem_base = nodes_elmnts(k+j*nsize)
             elem_target = nodes_elmnts(l+j*nsize)
             do n = 1, esize
                num_node = elmnts(esize*elem_base+n-1)
                do m = 0, nnodes_elmnts(num_node)-1
                   if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
                      connectivity = connectivity + 1
                   end if
                end do
             end do

             if ( connectivity >=  ncommonnodes) then

                is_neighbour = .false.

                do m = 0, xadj(nodes_elmnts(k+j*nsize))
                   if ( .not.is_neighbour ) then
                      if ( adjncy(nodes_elmnts(k+j*nsize)*sup_neighbour+m) == nodes_elmnts(l+j*nsize) ) then
                         is_neighbour = .true.

                      end if
                   end if
                end do
                if ( .not.is_neighbour ) then
                   adjncy(nodes_elmnts(k+j*nsize)*sup_neighbour+xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)
                   
                   xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
                   if (xadj(nodes_elmnts(k+j*nsize))>sup_neighbour) stop 'ERROR : too much neighbours per element, modify the mesh.'

                   adjncy(nodes_elmnts(l+j*nsize)*sup_neighbour+xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)

                   xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
                   if (xadj(nodes_elmnts(l+j*nsize))>sup_neighbour) stop 'ERROR : too much neighbours per element, modify the mesh.'
                end if
             end if
          end do
       end do
    end do

    max_neighbour = maxval(xadj)

    ! making adjacency arrays compact (to be used for partitioning)
    do i = 0, nelmnts-1
       k = xadj(i)
       xadj(i) = nb_edges
       do j = 0, k-1
          adjncy(nb_edges) = adjncy(i*sup_neighbour+j)
          nb_edges = nb_edges + 1
       end do
    end do

    xadj(nelmnts) = nb_edges


  end subroutine mesh2dual_ncommonnodes_fault

!-----------
!  subroutine write_fault_database_iface(IIN_database, iproc, nelmnts, &
!                                      glob2loc_elmnts, part)
!
!  integer, intent(in)  :: IIN_database
!  integer, intent(in)  :: iproc
!  integer(long), intent(in) :: nelmnts
!  integer, dimension(0:nelmnts-1), intent(in)  :: part,glob2loc_elmnts
!
!  integer, dimension(:), allocatable :: ispec1, ispec2, iface1, iface2 
!  integer  :: i,iflt,ispec_fault 
!  integer  :: nspec_fault_1,nspec_fault_2,nspec_fault
!  integer  :: k1,k2
!
!  do iflt=1,size(faults)
! 
!   ! check number of fault elements in this partition
!    nspec_fault_1 = count( part(faults(iflt)%ispec1-1) == iproc )
!    nspec_fault_2 = count( part(faults(iflt)%ispec2-1) == iproc )
!    print *, 'Fault # ',iflt
!    print *, '  ispec1 : ', nspec_fault_1
!    print *, '  ispec2 : ', nspec_fault_2
!    if (nspec_fault_1 /= nspec_fault_2) then
!      print *, 'Fatal error: Number of fault elements on ',iproc,' do not coincide. Abort.'
!      stop 
!    end if
!    nspec_fault = nspec_fault_1 
!
!    write(IIN_database,*) nspec_fault
!
!    if (nspec_fault==0) cycle 
!
!    allocate(ispec1(nspec_fault))
!    allocate(iface1(nspec_fault))
!    allocate(ispec2(nspec_fault))
!    allocate(iface2(nspec_fault))
!    k1 = 0
!    k2 = 0
!    do i = 1,faults(iflt)%nspec   
!      if ( part(faults(iflt)%ispec1(i)-1) == iproc ) then
!        k1 = k1 + 1
!        ispec1(k1)=glob2loc_elmnts(faults(iflt)%ispec1(i))
!        iface1(k1)=faults(iflt)%iface1(i)
!      endif
!      if ( part(faults(iflt)%ispec2(i)-1) == iproc ) then
!        k2 = k2 + 1
!        ispec2(k2)=glob2loc_elmnts(faults(iflt)%ispec2(i))
!        iface2(k2)=faults(iflt)%iface2(i)
!      endif
!    enddo
!
!   ! Writes ispec1 , ispec2 , iface1 , iface2
!    do i = 1,nspec_fault
!      write(IIN_database,*) ispec1(i), ispec2(i), iface1(i), iface2(i)
!    enddo
!   ! NOTE: the solver does not need ispec1 and ispec2 to be facing each other across the fault
!    deallocate(ispec1,ispec2,iface1,iface2)
!
!  enddo
!
!  end subroutine write_fault_database_iface



end module fault_scotch
