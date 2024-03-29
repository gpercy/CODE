!program pre_meshfem3D

module decompose_mesh_SCOTCH

  use part_decompose_mesh_SCOTCH
  use fault_scotch
  
  implicit none
  
  include './scotchf.h'

! number of partitions
  integer :: nparts ! e.g. 4 for partitioning for 4 CPUs or 4 processes

! mesh arrays
  integer(long) :: nspec
  integer, dimension(:,:), allocatable  :: elmnts
  integer, dimension(:,:), allocatable  :: mat
  integer, dimension(:), allocatable  :: part
  
  integer :: nnodes
  double precision, dimension(:,:), allocatable  :: nodes_coords
    
  integer, dimension(:), allocatable  :: xadj
  integer, dimension(:), allocatable  :: adjncy
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts
  integer, dimension(:), allocatable  :: elmnts_load

  integer, dimension(:), pointer  :: glob2loc_elmnts
  integer, dimension(:), pointer  :: glob2loc_nodes_nparts
  integer, dimension(:), pointer  :: glob2loc_nodes_parts
  integer, dimension(:), pointer  :: glob2loc_nodes

  integer, dimension(:), pointer  :: tab_size_interfaces, tab_interfaces
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces
  integer  ::  ninterfaces
  integer  :: my_ninterface
  
  integer(long)  :: nsize           ! Max number of elements that contain the same node.
  integer  :: nb_edges

  integer  :: ispec, inode
  integer  :: ngnod
  integer  :: max_neighbour         ! Real maximum number of neighbours per element
  integer(long)  :: sup_neighbour   ! Majoration of the maximum number of neighbours per element

  integer  :: ipart, nnodes_loc, nspec_loc
  integer  :: num_elmnt, num_node, num_mat

  ! boundaries
  integer  :: ispec2D
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top
  integer, dimension(:), allocatable :: ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin
  integer, dimension(:,:), allocatable :: nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top 

  ! moho surface (optional)
  integer :: nspec2D_moho
  integer, dimension(:), allocatable :: ibelm_moho
  integer, dimension(:,:), allocatable :: nodes_ibelm_moho
  
  character(len=256)  :: prname

  logical, dimension(:), allocatable :: mask_nodes_elmnts
  integer, dimension(:), allocatable :: used_nodes_elmnts

  double precision, dimension(SCOTCH_GRAPHDIM)  :: scotchgraph
  double precision, dimension(SCOTCH_STRATDIM)  :: scotchstrat
  character(len=256), parameter :: scotch_strategy='b{job=t,map=t,poli=S,sep=h{pass=30}}'
  integer  :: ierr,idummy
  
  !pll
  double precision , dimension(:,:), allocatable :: mat_prop
  integer :: count_def_mat,count_undef_mat,imat
  character (len=30), dimension(:,:), allocatable :: undef_mat_prop

! default mesh file directory
  character(len=256) :: localpath_name    
  character(len=256) :: outputpath_name 

  integer :: q_flag,aniso_flag,idomain_id
  double precision :: vp,vs,rho

  contains
  
  !----------------------------------------------------------------------------------------------
  ! reads in mesh files
  !----------------------------------------------------------------------------------------------
  subroutine read_mesh_files

  ! sets number of nodes per element
    ngnod = esize

  ! reads node coordinates
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file',&
          status='old', form='formatted', iostat = ierr)
    if( ierr /= 0 ) then
      print*,'could not open file:',localpath_name(1:len_trim(localpath_name))//'/nodes_coords_file'
      stop 'error file open'
    endif
    read(98,*) nnodes
    allocate(nodes_coords(3,nnodes))
    do inode = 1, nnodes
    ! format: #id_node #x_coordinate #y_coordinate #z_coordinate
      read(98,*) num_node, nodes_coords(1,num_node), nodes_coords(2,num_node), nodes_coords(3,num_node)
    !if(num_node /= inode)  stop "ERROR : Invalid nodes_coords file."
    end do
    close(98)
    print*, 'total number of nodes: '
    print*, '  nnodes = ', nnodes 

  ! reads mesh elements indexing 
  !(CUBIT calls this the connectivity, guess in the sense that it connects with the points index in 
  ! the global coordinate file "nodes_coords_file"; it doesn't tell you which point is connected with others)
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/mesh_file', &
          status='old', form='formatted')
    read(98,*) nspec
    allocate(elmnts(esize,nspec))
    do ispec = 1, nspec
      ! format: # element_id  #id_node1 ... #id_node8

      ! note: be aware that here we can have different node ordering for a cube element;
      !          the ordering from Cubit files might not be consistent for multiple volumes, or uneven, unstructured grids
      !         
      !          guess here it assumes that spectral elements ordering is like first at the bottom of the element, anticlock-wise, i.e. 
      !             point 1 = (0,0,0), point 2 = (0,1,0), point 3 = (1,1,0), point 4 = (1,0,0)
      !          then top (positive z-direction) of element 
      !             point 5 = (0,0,1), point 6 = (0,1,1), point 7 = (1,1,1), point 8 = (1,0,1)

      !read(98,*) num_elmnt, elmnts(5,num_elmnt), elmnts(1,num_elmnt),elmnts(4,num_elmnt), elmnts(8,num_elmnt), &
      !      elmnts(6,num_elmnt), elmnts(2,num_elmnt), elmnts(3,num_elmnt), elmnts(7,num_elmnt)

      read(98,*) num_elmnt, elmnts(1,num_elmnt), elmnts(2,num_elmnt),elmnts(3,num_elmnt), elmnts(4,num_elmnt), &
            elmnts(5,num_elmnt), elmnts(6,num_elmnt), elmnts(7,num_elmnt), elmnts(8,num_elmnt)

      if((num_elmnt > nspec) .or. (num_elmnt < 1) )  stop "ERROR : Invalid mesh file."

        
      !outputs info for each element to see ordering
      !print*,'ispec: ',ispec
      !print*,'  ',num_elmnt, elmnts(5,num_elmnt), elmnts(1,num_elmnt),elmnts(4,num_elmnt), elmnts(8,num_elmnt), &
      !      elmnts(6,num_elmnt), elmnts(2,num_elmnt), elmnts(3,num_elmnt), elmnts(7,num_elmnt)    
      !print*,'elem:',num_elmnt
      !do i=1,8
      !  print*,' i ',i,'val :',elmnts(i,num_elmnt),&
      !    nodes_coords(1,elmnts(i,num_elmnt)),nodes_coords(2,elmnts(i,num_elmnt)),nodes_coords(3,elmnts(i,num_elmnt))
      !enddo
      !print*
          
    end do
    close(98)
    print*, 'total number of spectral elements:'
    print*, '  nspec = ', nspec

  ! reads material associations
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/materials_file', &
          status='old', form='formatted')
    allocate(mat(2,nspec))
    do ispec = 1, nspec
      ! format: # id_element #flag
      ! note: be aware that elements may not be sorted in materials_file
      read(98,*) num_mat,mat(1,num_mat) !mat(1,ispec)!, mat(2,ispec) 
      if((num_mat > nspec) .or. (num_mat < 1) ) stop "ERROR : Invalid mat file."
    end do
    close(98)

  ! TODO:
  ! must be changed, if  mat(1,i) < 0  1 == interface , 2 == tomography
    mat(2,:) = 1
    
  ! reads material definitions
  !
  ! note: format of nummaterial_velocity_file must be
  !
  ! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_flag  #(7)anisotropy_flag
  !
  ! where
  !     material_domain_id : 1=acoustic / 2=elastic / 3=poroelastic
  !     material_id               : number of material/volume
  !     rho                           : density
  !     vp                             : P-velocity
  !     vs                             : S-velocity
  !     Q_flag                      : 0=no attenuation/1=IATTENUATION_SEDIMENTS_40, 2=..., 13=IATTENUATION_BEDROCK
  !     anisotropy_flag        : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
    count_def_mat = 0
    count_undef_mat = 0
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_velocity_file',&
          status='old', form='formatted')
    ! note: format #material_domain_id #material_id #...      
    read(98,*,iostat=ierr) idummy,num_mat
    print *,'materials:'
    ! counts materials (defined/undefined)
    do while (ierr == 0)
       print*, '  num_mat = ',num_mat
       if(num_mat /= -1) then 
          count_def_mat = count_def_mat + 1        
       else
          count_undef_mat = count_undef_mat + 1
       end if
       read(98,*,iostat=ierr) idummy,num_mat
    end do
    close(98)
    print*, '  defined = ',count_def_mat, 'undefined = ',count_undef_mat
    ! check with material flags
    if( count_def_mat > 0 .and. maxval(mat(1,:)) > count_def_mat ) then
      print*,'error material definitions:'
      print*,'  materials associated in materials_file:',maxval(mat(1,:))
      print*,'  bigger than defined materials in nummaterial_velocity_file:',count_def_mat
      stop 'error materials'
    endif
    allocate(mat_prop(6,count_def_mat))
    allocate(undef_mat_prop(6,count_undef_mat))
    ! reads in defined material properties
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/nummaterial_velocity_file', &
          status='old', form='formatted')
    do imat=1,count_def_mat
       ! material definitions
       !
       ! format: note that we save the arguments in a slightly different order in mat_prop(:,:)
       !              #(6) material_domain_id #(0) material_id  #(1) rho #(2) vp #(3) vs #(4) Q_flag #(5) anisotropy_flag
       !
       read(98,*) idomain_id,num_mat,rho,vp,vs,q_flag,aniso_flag
       !read(98,*) num_mat, mat_prop(1,num_mat),mat_prop(2,num_mat),&
       !           mat_prop(3,num_mat),mat_prop(4,num_mat),mat_prop(5,num_mat)
       mat_prop(1,num_mat) = rho
       mat_prop(2,num_mat) = vp
       mat_prop(3,num_mat) = vs
       mat_prop(4,num_mat) = q_flag
       mat_prop(5,num_mat) = aniso_flag
       mat_prop(6,num_mat) = idomain_id
       
       if(num_mat < 0 .or. num_mat > count_def_mat)  stop "ERROR : Invalid nummaterial_velocity_file file."    

       !checks attenuation flag with integer range as defined in constants.h like IATTENUATION_SEDIMENTS_40, ....
       if( int(mat_prop(4,num_mat)) > 13 ) then
          stop 'wrong attenuation flag in mesh: too large, not supported yet - check with constants.h'
       endif
    end do
    ! reads in undefined material properties
    do imat=1,count_undef_mat
       read(98,'(6A30)') undef_mat_prop(6,imat),undef_mat_prop(1,imat),undef_mat_prop(2,imat),&
                        undef_mat_prop(3,imat),undef_mat_prop(4,imat),undef_mat_prop(5,imat)
    end do
    close(98)

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmin', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_xmin = 0
    else
      read(98,*) nspec2D_xmin
    endif
    allocate(ibelm_xmin(nspec2D_xmin))
    allocate(nodes_ibelm_xmin(4,nspec2D_xmin))
    do ispec2D = 1,nspec2D_xmin 
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      ! note: ordering for CUBIT seems such that the normal of the face points outward of the element the face belongs to;
      !         in other words, nodes are in increasing order such that when looking from within the element outwards, 
      !         they are ordered clockwise
      !
      !          doesn't necessarily have to start on top-rear, then bottom-rear, bottom-front, and finally top-front i.e.: 
      !          point 1 = (0,1,1), point 2 = (0,1,0), point 3 = (0,0,0), point 4 = (0,0,1)
      read(98,*) ibelm_xmin(ispec2D), nodes_ibelm_xmin(1,ispec2D), nodes_ibelm_xmin(2,ispec2D), &
            nodes_ibelm_xmin(3,ispec2D), nodes_ibelm_xmin(4,ispec2D)

      !outputs info for each element for check of ordering          
      !print*,'ispec2d:',ispec2d
      !print*,'  xmin:', ibelm_xmin(ispec2D), nodes_ibelm_xmin(1,ispec2D), nodes_ibelm_xmin(2,ispec2D), &
      !      nodes_ibelm_xmin(3,ispec2D), nodes_ibelm_xmin(4,ispec2D)     
      !do i=1,4
      !  print*,'i',i,'val:',ibelm_xmin(ispec2d),nodes_coords(1,nodes_ibelm_xmin(i,ispec2D)), &
      !      nodes_coords(2,nodes_ibelm_xmin(i,ispec2D)),nodes_coords(3,nodes_ibelm_xmin(i,ispec2D))
      !enddo
      !print*
    end do
    close(98)
    print*, 'absorbing boundaries:'
    print*, '  nspec2D_xmin = ', nspec2D_xmin 

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_xmax', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_xmax = 0
    else
      read(98,*) nspec2D_xmax
    endif
    allocate(ibelm_xmax(nspec2D_xmax))
    allocate(nodes_ibelm_xmax(4,nspec2D_xmax))
    do ispec2D = 1,nspec2D_xmax
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_xmax(ispec2D), nodes_ibelm_xmax(1,ispec2D), nodes_ibelm_xmax(2,ispec2D), &
            nodes_ibelm_xmax(3,ispec2D), nodes_ibelm_xmax(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_xmax = ', nspec2D_xmax

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymin', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_ymin = 0
    else
      read(98,*) nspec2D_ymin
    endif
    allocate(ibelm_ymin(nspec2D_ymin))
    allocate(nodes_ibelm_ymin(4,nspec2D_ymin))
    do ispec2D = 1,nspec2D_ymin 
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face   
      read(98,*) ibelm_ymin(ispec2D), nodes_ibelm_ymin(1,ispec2D), nodes_ibelm_ymin(2,ispec2D),  &
            nodes_ibelm_ymin(3,ispec2D), nodes_ibelm_ymin(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_ymin = ', nspec2D_ymin 

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_ymax', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_ymax = 0
    else
      read(98,*) nspec2D_ymax
    endif
    allocate(ibelm_ymax(nspec2D_ymax))
    allocate(nodes_ibelm_ymax(4,nspec2D_ymax))
    do ispec2D = 1,nspec2D_ymax 
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face  
      read(98,*) ibelm_ymax(ispec2D), nodes_ibelm_ymax(1,ispec2D), nodes_ibelm_ymax(2,ispec2D),  &
            nodes_ibelm_ymax(3,ispec2D), nodes_ibelm_ymax(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_ymax = ', nspec2D_ymax

  ! reads in absorbing boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/absorbing_surface_file_bottom', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_bottom = 0
    else
      read(98,*) nspec2D_bottom
    endif
    allocate(ibelm_bottom(nspec2D_bottom))
    allocate(nodes_ibelm_bottom(4,nspec2D_bottom))
    do ispec2D = 1,nspec2D_bottom 
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face   
      read(98,*) ibelm_bottom(ispec2D), nodes_ibelm_bottom(1,ispec2D), nodes_ibelm_bottom(2,ispec2D), &
            nodes_ibelm_bottom(3,ispec2D), nodes_ibelm_bottom(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_bottom = ', nspec2D_bottom 

  ! reads in free_surface boundary files
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/free_surface_file', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_top = 0
    else
      read(98,*) nspec2D_top
    endif
    allocate(ibelm_top(nspec2D_top))
    allocate(nodes_ibelm_top(4,nspec2D_top))
    do ispec2D = 1,nspec2D_top 
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_top(ispec2D), nodes_ibelm_top(1,ispec2D), nodes_ibelm_top(2,ispec2D), &
             nodes_ibelm_top(3,ispec2D), nodes_ibelm_top(4,ispec2D)
    end do
    close(98)
    print*, '  nspec2D_top = ', nspec2D_top

  ! reads in moho_surface boundary files (optional)
    open(unit=98, file=localpath_name(1:len_trim(localpath_name))//'/moho_surface_file', &
          status='old', form='formatted',iostat=ierr)
    if( ierr /= 0 ) then
      nspec2D_moho = 0
    else
      read(98,*) nspec2D_moho
    endif
    allocate(ibelm_moho(nspec2D_moho))
    allocate(nodes_ibelm_moho(4,nspec2D_moho))
    do ispec2D = 1,nspec2D_moho
      ! format: #id_(element containing the face) #id_node1_face .. #id_node4_face
      read(98,*) ibelm_moho(ispec2D), nodes_ibelm_moho(1,ispec2D), nodes_ibelm_moho(2,ispec2D), &
             nodes_ibelm_moho(3,ispec2D), nodes_ibelm_moho(4,ispec2D)
    end do
    close(98)
    if( nspec2D_moho > 0 ) print*, '  nspec2D_moho = ', nspec2D_moho
  print*, 'test2'
    call read_fault_files(localpath_name)
  print*, 'test3'
    if (allocated(faults)) then
      print*, 'test4' 
      call save_nodes_coords(nodes_coords,nnodes)
      print*, 'test5'
      call close_faults(nodes_coords,elmnts,nspec,nnodes,esize)    
      print*, 'test6'
    end if 

  print*, 'test1'
  end subroutine read_mesh_files
  
  !----------------------------------------------------------------------------------------------
  ! checks valence of nodes
  !----------------------------------------------------------------------------------------------
  
  subroutine check_valence
  
    allocate(mask_nodes_elmnts(nnodes))
    allocate(used_nodes_elmnts(nnodes))
    mask_nodes_elmnts(:) = .false.
    used_nodes_elmnts(:) = 0
    do ispec = 1, nspec
      do inode = 1, ESIZE
        mask_nodes_elmnts(elmnts(inode,ispec)) = .true.
        used_nodes_elmnts(elmnts(inode,ispec)) = used_nodes_elmnts(elmnts(inode,ispec)) + 1
      enddo
    enddo
    print *, 'nodes valence: '
    print *, '  min = ',minval(used_nodes_elmnts(:)),'max = ', maxval(used_nodes_elmnts(:))
    do inode = 1, nnodes
      if (.not. mask_nodes_elmnts(inode)) then
        stop 'ERROR : nodes not used.'
      endif
    enddo
    nsize = maxval(used_nodes_elmnts(:))
    sup_neighbour = ngnod * nsize - (ngnod + (ngnod/2 - 1)*nfaces)
    print*, '  nsize = ',nsize, 'sup_neighbour = ', sup_neighbour

  end subroutine check_valence

  !----------------------------------------------------------------------------------------------
  ! divides model into partitions using scotch library functions
  !----------------------------------------------------------------------------------------------
  
  subroutine scotch_partitioning

    implicit none

    elmnts(:,:) = elmnts(:,:) - 1

    ! determines maximum neighbors based on 1 common node
    allocate(xadj(1:nspec+1))
    allocate(adjncy(1:sup_neighbour*nspec))
    allocate(nnodes_elmnts(1:nnodes))
    allocate(nodes_elmnts(1:nsize*nnodes))    
    call mesh2dual_ncommonnodes(nspec, nnodes, nsize, sup_neighbour, elmnts, xadj, adjncy, nnodes_elmnts, &
         nodes_elmnts, max_neighbour, 1)
    print*, 'mesh2dual: '
    print*, '  max_neighbour = ',max_neighbour

    nb_edges = xadj(nspec+1)

  ! allocates & initializes partioning of elements
    allocate(part(1:nspec))
    part(:) = -1

  ! initializes
  ! elements load array
    allocate(elmnts_load(1:nspec))
    
    ! uniform load
    elmnts_load(:) = 1 
    
    ! in case of acoustic/elastic simulation, weights elements accordingly
    call acoustic_elastic_load(elmnts_load,nspec,count_def_mat,mat(1,:),mat_prop)
    
  ! SCOTCH partitioning
    call scotchfstratinit (scotchstrat(1), ierr)
     if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot initialize strat'
    endif

    call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), ierr)
     if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot build strat'
    endif

    call scotchfgraphinit (scotchgraph (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot initialize graph'
    endif

    ! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
    ! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
    !                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
    !                    #(6) vertex_load_array (optional) #(7) vertex_label_array
    !                    #(7) number_of_arcs                    #(8) adjacency_array 
    !                    #(9) arc_load_array (optional)      #(10) ierror
    call scotchfgraphbuild (scotchgraph (1), 0, nspec, &
                          xadj (1), xadj (1), &
                          elmnts_load (1), xadj (1), &
                          nb_edges, adjncy (1), &
                          adjncy (1), ierr)

    ! w/out element load, but adjacency array
    !call scotchfgraphbuild (scotchgraph (1), 0, nspec, &
    !                      xadj (1), xadj (1), &
    !                      xadj (1), xadj (1), &
    !                      nb_edges, adjncy (1), &
    !                      adjncy (1), ierr)
                          
                          
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot build graph'
    endif

    call scotchfgraphcheck (scotchgraph (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Invalid check'
    endif

    call scotchfgraphpart (scotchgraph (1), nparts, scotchstrat(1),part(1),ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot part graph'
    endif

    call scotchfgraphexit (scotchgraph (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot destroy graph'
    endif

    call scotchfstratexit (scotchstrat(1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot destroy strat'
    endif


  ! re-partitioning puts acoustic-elastic coupled elements into same partition
  !  integer  :: nfaces_coupled
  !  integer, dimension(:,:), pointer  :: faces_coupled
  !    call acoustic_elastic_repartitioning (nspec, nnodes, elmnts, &
  !                   count_def_mat, mat(1,:) , mat_prop, &
  !                   sup_neighbour, nsize, &
  !                   nparts, part, nfaces_coupled, faces_coupled)
 
  ! move all fault elements to the same partition (proc=0)
    call fault_repartition (nspec, nnodes, elmnts, nsize, nparts, part, esize)

     
  ! re-partitioning puts moho-surface coupled elements into same partition
    call moho_surface_repartitioning (nspec, nnodes, elmnts, &
                     sup_neighbour, nsize, nparts, part, &
                     nspec2D_moho,ibelm_moho,nodes_ibelm_moho )
   

  ! local number of each element for each partition
    call Construct_glob2loc_elmnts(nspec, part, glob2loc_elmnts,nparts)


  ! local number of each node for each partition
    call Construct_glob2loc_nodes(nspec, nnodes,nsize, nnodes_elmnts, nodes_elmnts, part, &
         glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, nparts)


  ! mpi interfaces 
    ! acoustic/elastic boundaries WILL BE SEPARATED into different MPI partitions
    call Construct_interfaces(nspec, sup_neighbour, part, elmnts, &
                             xadj, adjncy, tab_interfaces, &
                             tab_size_interfaces, ninterfaces, &
                             nparts)

    !or: uncomment if you want acoustic/elastic boundaries NOT to be separated into different MPI partitions
    !call Construct_interfaces_no_ac_el_sep(nspec, sup_neighbour, part, elmnts, &
    !                          xadj, adjncy, tab_interfaces, &
    !                          tab_size_interfaces, ninterfaces, &
    !                          count_def_mat, mat_prop(3,:), mat(1,:), nparts)


  end subroutine scotch_partitioning

 
  !----------------------------------------------------------------------------------------------
  ! writes out new Databases files for each partition
  !----------------------------------------------------------------------------------------------
  
  subroutine write_mesh_databases

    allocate(my_interfaces(0:ninterfaces-1))
    allocate(my_nb_interfaces(0:ninterfaces-1))

    ! writes out Database file for each partition

    do ipart = 0, nparts-1

       ! opens output file
       write(prname, "(i6.6,'_Database')") ipart
       open(unit=15,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname,&
            status='unknown', action='write', form='formatted', iostat = ierr)
       if( ierr /= 0 ) then
        print*,'error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
        print*
        print*,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
        stop 'error file open Database'
       endif

       ! gets number of nodes 
       call write_glob2loc_nodes_database(15, ipart, nnodes_loc, nodes_coords, &
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 1) 

       ! gets number of spectral elements                           
       call write_partition_database(15, ipart, nspec_loc, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, ngnod, 1)

       ! writes out node coordinate locations 
       write(15,*) nnodes_loc
       
       call write_glob2loc_nodes_database(15, ipart, nnodes_loc, nodes_coords,&
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 2)

       call write_material_properties_database(15,count_def_mat,count_undef_mat, &
                                  mat_prop, undef_mat_prop) 
        
       ! writes out spectral element indices 
       write(15,*) nspec_loc
       
       call write_partition_database(15, ipart, nspec_loc, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, ngnod, 2)

       ! writes out absorbing/free-surface boundaries
       call write_boundaries_database(15, ipart, nspec, nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                  nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                  ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                  ibelm_ymax, ibelm_bottom, ibelm_top, &
                                  nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                  nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part)

       ! gets number of MPI interfaces                           
       call Write_interfaces_database(15, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                  my_ninterface, my_interfaces, my_nb_interfaces, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, 1, nparts)

       ! writes out MPI interfaces elements
       write(15,*) my_ninterface, maxval(my_nb_interfaces)
       
       call Write_interfaces_database(15, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                  my_ninterface, my_interfaces, my_nb_interfaces, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, 2, nparts)       

       ! writes out moho surface (optional) 
       call write_moho_surface_database(15, ipart, nspec, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, &
                                  nspec2D_moho,ibelm_moho,nodes_ibelm_moho)
        
       close(15)
      
       if (.not. allocated(faults)) cycle 
          ! write fault database
          write(prname, "(i6.6,'_Database_fault')") ipart
          open(unit=16,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname,&
               status='unknown', action='write', form='formatted', iostat = ierr)
          if( ierr /= 0 ) then
            print*,'error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
            print*
            print*,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
            stop 
          endif
          call write_fault_database(16, ipart, nspec, &
                                    glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                    glob2loc_nodes, part)
!          close(16)
         
         ! write nodes coordinates with fault open crack  
!          write(prname, "(i6.6,'_Database_fault_nodes')") ipart
!          open(unit=16,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname,&
!               status='unknown', action='write', form='formatted', iostat = ierr)
!          if( ierr /= 0 ) then
!            print*,'error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
!            print*
!            print*,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
!            stop 
!          endif
   
          write(16,*) nnodes_loc
          call write_glob2loc_nodes_database(16, ipart, nnodes_loc, nodes_coords_open,&
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, nnodes, 2)

          close(16)

    end do
    print*, 'partitions: '
    print*, '  num = ',nparts
    print*
    print*, 'Databases files in directory: ',outputpath_name(1:len_trim(outputpath_name))
    print*, 'finished successfully'
    print*

  end subroutine write_mesh_databases
  
!end program pre_meshfem3D

end module decompose_mesh_SCOTCH

