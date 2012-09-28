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

  public :: read_fault_files, fault_repartition, close_faults, write_fault_database, &
            save_nodes_coords, nodes_coords_open, faults

CONTAINS 
!==========================================================================================

  Subroutine read_fault_files(localpath_name)

  character(len=256),intent(in) :: localpath_name    
  integer :: nbfaults, iflt, ier 

  open(101,file='../DATA/FAULT/Par_file_faults',status='old',action='read',iostat=ier)
  if (ier==0) then 
    read(101,*) nbfaults
  else
    nbfaults = 0
    print *, 'Par_file_faults not found: assume no faults'
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
  enddo
  do e=1,f%nspec
    read(101,*) f%ispec2(e),f%inodes2(:,e)
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

  Subroutine fault_repartition (nelmnts, nnodes, elmnts, nsize, &
                        nproc, part, esize)

!     part : iproc number of processor partioned. It will altered patching fault elements into the same partion.  
!     Part, once is altered , will be input for write_partition_database.

!INPUTS
  integer(long), intent(in) :: nelmnts,nsize
  integer, intent(in)  :: nnodes, nproc, esize 
  integer, dimension(0:esize*nelmnts-1), intent(in) :: elmnts
!OUTPUTS :
  integer, dimension(0:nelmnts-1), intent(inout)    :: part

  integer, dimension(0:nnodes-1)                 :: nnodes_elmnts
  integer, dimension(0:nsize*nnodes-1)           :: nodes_elmnts
  integer  :: max_neighbour       

  integer  :: i,j, ipart,nproc_null,nproc_null_final
  integer  :: k1, k2, k,e,iflt,inode
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

end module fault_scotch
