module IOLib
   use MeshLib
   implicit none
   type :: IOConfig
      ! DON'T CHANGE WITHOUT KNOWING WHAT YOU'RE DOING
      ! OTHER SUBROUTINES HERE USE THESE DIRECTORIES IN CONCATENATING THEIR FILENAMES
      ! output directory, where all output files are dumped
      character(len=9) :: outputDir = "./output/"
      ! input director
      character(len=8) :: inputDir = "./input/"
   end type IOConfig
contains

   subroutine generateLevelMeshElementsFile(base_mesh, meshLevel)
      implicit none
      type(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel
      type(IOConfig) :: io
      character(len=128) :: elemsFname
      integer :: i, j
      type(Tetrahedron) :: tet
      integer :: vertices(4)
      !   type(LevelMesh) :: mesh

      ! output file to save elements per level
      if(meshLevel<10) then
         write(elemsFname, "(A9,A9,I1,A4)") io%outputDir, "elements_", meshLevel, ".dat"
      elseif(meshLevel<100) then
         write(elemsFname, "(A9,A9,I2,A4)") io%outputDir, "elements_", meshLevel, ".dat"
      elseif(meshLevel<1000) then
         write(elemsFname, "(A9,A9,I3,A4)") io%outputDir, "elements_", meshLevel, ".dat"
      elseif(meshLevel<10000) then
         write(elemsFname, "(A9,A9,I4,A4)") io%outputDir, "elements_", meshLevel, ".dat"
      else
         write(*,*) 'ERROR. in generating mesh elements file. Can not do more than 9999 output files.'
         stop
      endif

      ! For any level, you want to include the elements from all
      ! the levels before it
      open(unit=19, file=elemsFname, access="sequential")
      ! loop from meshLevel 1 to meshLevel
      ! access the elements and print them in order. the globalIds should be in asc
      ! order, incrementing by 1 with each step

      do i = 1, meshLevel
         ! over the elements in this level
         do j = 1, base_mesh%levelMeshes(i)%numElements
            ! extract tet
            ! to get global Id of tet at this position in order to get tet from base_mesh
            tet = base_mesh%tetrahedrons(base_mesh%levelMeshes(i)%elements(j))
            ! Format: elemGlobalId, level, levelId, nodes' connectivity
            vertices = (/tet%vertices(1)%globalId, tet%vertices(2)%globalId,&
               tet%vertices(3)%globalId, tet%vertices(4)%globalId/)
            write(19, *) tet%globalId, tet%level, tet%elementLevelId, vertices
         end do
      end do
      close(19)
   end subroutine generateLevelMeshElementsFile

   subroutine generateNodesFile(base_mesh)
      implicit none
      type(BaseMesh), intent(in) :: base_mesh
      type(IOConfig) :: io
      character(len=128) :: nodesFname
      integer :: i
      type(Vertex) :: vtx

      ! output file to save all nodes of base_mesh
      nodesFname = io%outputDir // "nodes.dat"
      open(unit=19, file=nodesFname, access="sequential")

      ! over all nodes in basemesh
      do i = 1, base_mesh%numVertices
         ! format: globalNodeId, parentLevel, x, y, z
         vtx = base_mesh%v2pe(i)
         write(19, *) vtx%globalId, vtx%elementLevel, vtx%x, vtx%y, vtx%z
      end do
      close(19)
   end subroutine generateNodesFile

   subroutine generateLevelMeshHangingNodesFile(base_mesh, meshLevel)
      implicit none
      type(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel
      type(Vertex) :: vtx
      type(IOConfig) :: io
      integer :: i
      type(Vertex), allocatable :: hangingNodes(:)
      character(len=128) :: nodesFname

      hangingNodes = base_mesh%getLevelMeshHangingNodes(meshLevel)

      ! output file to save all nodes of base_mesh
      if(meshLevel<10) then
         write(nodesFname, "(A9,A13,I1,A4)") io%outputDir, "hangingnodes_", meshLevel, ".dat"
      elseif(meshLevel<100) then
         write(nodesFname, "(A9,A13,I2,A4)") io%outputDir, "hangingnodes_", meshLevel, ".dat"
      elseif(meshLevel<1000) then
         write(nodesFname, "(A9,A13,I3,A4)") io%outputDir, "hangingnodes_", meshLevel, ".dat"
      elseif(meshLevel<10000) then
         write(nodesFname, "(A9,A13,I4,A4)") io%outputDir, "hangingnodes_", meshLevel, ".dat"
      else
         write(*,*) 'ERROR. in generating hangingnodes file. Can not do more than 9999 output files.'
         stop
      endif

      open(unit=19, file=nodesFname, access="sequential")

      ! over all nodes in returned hanging Nodes array
      do i = 1, size(hangingNodes)
         ! format: globalNodeId, parentLevel, x, y, z
         vtx = hangingNodes(i)
         write(19, *) vtx%globalId, vtx%elementLevel, vtx%x, vtx%y, vtx%z
      end do
      close(19)
   end subroutine generateLevelMeshHangingNodesFile

   subroutine generateNodesInterpolationFile(base_mesh)
      implicit none
      type(BaseMesh), intent(in) :: base_mesh
      type(Vertex) :: vtx
      integer :: numLevelNodes, i, meshLevel
      character(len=128) :: intpFname
      type(IOConfig) :: io

      intpFname = io%outputDir // "intpdata.dat"

      ! structure of interpolation data file
      ! for each level:
      ! - level's numNodes
      ! for each node on level
      ! -global node iD
      ! -number of parents (0 for coarsest mesh i.e. level 1)
      ! - parent global node IDs in ASC order (not written if no parents)
      !(Each node of finer meshes (meshLevel>1) that also exists on the coarser mesh now has
      ! one parent - itself with interpolation weight 1.0)
      ! - associated interpolation weights (we use a linear interpolation, getting the
      ! midpoints of parentEdges i.e. v_child = 0.5v1 + 0.5v2)
      open(unit=73,file=intpFname,access='sequential')

      ! over the levels
      do meshLevel = 1, base_mesh%numLevels
         numLevelNodes = base_mesh%getLevelMeshNumNodes(meshLevel)
         write(73, *) numLevelNodes
         ! over the level nodes
         do i = 1, numLevelNodes
            vtx = base_mesh%v2pe(i)
            write(73, *) vtx%globalId
            ! for coarsest mesh
            if ( meshLevel == 1 ) then
               ! zero parent
               write(73,*) 0
            end if

            if (meshLevel > 1) then
               ! finer meshes
               ! number of parent & parent list
               ! for nodes from any prior coarse mesh that appear in finer meshes - 1 parent
               if ( vtx%elementLevel < (meshLevel-1) ) then
                  write(73,*) 1
                  ! parent - itself
                  write(73, *) vtx%globalId
                  ! associated interpolation weight
                  write(73,*) 1.0
                  ! ACTUAL new nodes from (meshLevel-1) to form elements on `meshLevel`
                  ! ACTUAL nodes that make `meshLevel` finer
               else if (vtx%elementLevel == (meshLevel-1)) then
                  ! num of parents
                  write(73,*) 2
                  ! parent list
                  write(73,*) vtx%parentEdge%vertices(1), vtx%parentEdge%vertices(2)
                  ! weights
                  write(73,*) 0.5, 0.5
               end if
            end if
         end do
      end do

   end subroutine generateNodesInterpolationFile

   subroutine generateLevelMeshViewFile(base_mesh, meshLevel)
      implicit none
      type(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel
      type(IOConfig) :: io
      character(len=128) :: viewFname
      integer :: numElements, numVertices, i, j
      type(Vertex) :: vtx
      !   type(Vertex), allocatable :: levelVertices(:)
      viewFname = io%outputDir // 'view_mesh_'!Base name of output file. Must have length 4

      if(meshLevel<10) then
         write(viewFname, "(A19,I1,A4)") viewFname, meshLevel,".vtu"
      elseif(meshLevel<100) then
         write(viewFname, "(A19,I2,A4)") viewFname, meshLevel,".vtu"
      elseif(meshLevel<1000) then
         write(viewFname, "(A19,I3,A4)") viewFname, meshLevel,".vtu"
      elseif(meshLevel<10000) then
         write(viewFname, "(A19,I4,A4)") viewFname, meshLevel,".vtu"
      else
         write(*,*) 'ERROR. In generating mesh view file. Can not do more than 9999 output files.'
         stop
      endif

      open(unit=99,file=viewFname,access='sequential')
      write(*,*) 'Saving PARAVIEW data to ',viewFname

      ! Write out top header:
      write(99,1000)

      ! The number of elements and nodes in level `meshLevel` should be the sum of
      ! the elements from 1 to meshLevel. Similarly, nodes that form elements in meshLevel
      ! are L(meshLevel-1) nodes or from a higher level => we take the node of elementLevel
      ! from 0 to meshLevel -1. to get numElements, we loop sum over the levels and add
      ! each level mesh's numElements
      numElements = 0
      do i = 1, meshLevel
         numElements = numElements + base_mesh%levelMeshes(i)%numElements
      end do

      ! for numVertices, we have to count the vertices having elementLevels from 0 to meshLevel-1
      numVertices = 0
      do i = 1, base_mesh%numVertices
         vtx = base_mesh%v2pe(i)
         if ( vtx%elementLevel < meshLevel ) then
            numVertices = numVertices + 1
         end if
         if(vtx%elementLevel >= meshLevel) exit
      end do

      write(99,1010) numVertices, numElements  ! Start Mesh/Piece Section
      !numnp1, numnp2: node numbers -> replace numnp1+numnp2 by your number of nodes
      !numel(1): number of elements

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(99,1020)                           ! Start Point/Node data
      write(99,1030) 'Float64','nodes',3
      ! loop over base_mesh%v2pe for numVertices times
      do i = 1,numVertices
         vtx = base_mesh%v2pe(i)
         write(99,5000) vtx%x
         write(99,5000) vtx%y
         write(99,5000) vtx%z
      end do
      write(99,1031)                           ! Close Node data
      write(99,1021)                           ! Close Points section
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



      !=============================================================
      write(99,1050)                           ! Start Cell Section

      write(99,1030) 'Int32','connectivity',1  ! Start Elements
      ! over base_mesh's tets for numElements times
      do i = 1, numElements
         ! Paraview's count starts at zero and so the -1 to offset the start of 1
         write(99,6000) base_mesh%tetrahedrons(i)%vertices(1)%globalId -1
         write(99,6000) base_mesh%tetrahedrons(i)%vertices(2)%globalId -1
         write(99,6000) base_mesh%tetrahedrons(i)%vertices(3)%globalId -1
         write(99,6000) base_mesh%tetrahedrons(i)%vertices(4)%globalId -1
      end do
      write(99,1031)                          ! Close Elements

      write(99,1030) 'Int32','offsets',1      ! Start Offsets
      !   j=nodes_per_element
      j=4
      do i = 1, numElements
         write(99,6000) j
         !  j = j + nodes_per_element
         j = j+4
      end do
      write(99,1031)                           ! Close Offsets

      ! Paraview's VTK files use integers to represent different element types
      write(99,1030) 'UInt8','types',1         ! Start Element types
      ! 4 node tetrahedron. Tet Library, no one expects a brick element to be sent
      do i = 1, numElements
         write(99,6000) 10
      end do
      write(99,1031)                             ! Close Element types

      write(99,1051)                             ! Close Cell Section
      !=============================================================

      write(99,1011)                             ! Close Mesh/Piece


      ! Close the XML file:
      write(99,1001)

1000  format('<?xml version="1.0"?>'/'<VTKFile type="UnstructuredGrid" version="0.1">'/'<UnstructuredGrid>')
1001  format('</UnstructuredGrid> </VTKFile>')


1010  format('<Piece NumberOfPoints="',i10,'" NumberOfCells="',i10,'">')
1011  format('</Piece>')


1020  format('<Points>')
1021  format('</Points>')

1030  format('<DataArray type="',a,'" Name="',a,'" NumberOfComponents="',i2,'" format="ascii">')
1031  format('</DataArray>')

1050  format('<Cells>')
1051  format('</Cells>')

5000  format(e15.5,' ',$)
6000  format(i6,' ',$)
   end subroutine generateLevelMeshViewFile
end module IOLib
