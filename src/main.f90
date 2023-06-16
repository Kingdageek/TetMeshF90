program main
   use MeshLib
   use helpers
   use IOLib
   ! use ftlHashSetIntModule
   implicit none
   type(BaseMesh) :: base_mesh
   type(Vertex) :: nodes(5), elementNodes(4)
   type(Tetrahedron), allocatable :: tet(:)
   type(MeshEdge) :: edge
   integer :: numNodes, currLevel, numElements, numLevelElements, i,j, domainId, vtxLevel
   type(Tetrahedron) :: neighbours(4)
   integer, allocatable :: intArr(:)
   type(Vertex), allocatable :: hangingNodes(:)
   integer :: nextLevel
   double precision :: vertices(5,3)
   double precision :: x, y, z
   ! type(ftlHashSetInt) :: intSet, intSet2, resultSet
   ! type(ftlHashSetIntIterator) :: intIter, resIter

   vertices(1,:) = (/0.0, 0.0, 0.0/)
   vertices(2,:) = (/1.0, -1.73, 0.0/)
   vertices(3,:) = (/2.0, 0.0, 0.0/)
   vertices(4,:) = (/1.0, -0.58, 1.63/)
   vertices(5,:) = (/2.67, -1.54, 1.09/)
   ! The entire domain being meshed is taken to be on Level 0
   domainId = 0

   ! To try to model a mesh of two tetrahedrons sharing a face
   allocate(tet(2))
   ! we need a BaseMesh struct. initialize the struct

   write(*,*) "initializing BaseMesh..."
   call base_mesh%initialize()
   write(*,*) "initializing LevelMesh..."
   ! initialize a levelMesh for level 1
   call base_mesh%initializeLevelMesh()

   write(*,*) "successfully initialized BaseMesh and LevelMesh..."

   write(*,*) "Creating nodes for 1st element..."

   ! For tracking
   currLevel = base_mesh%numLevels
   numNodes = base_mesh%numVertices
   numElements = base_mesh%numElements
   numLevelElements = base_mesh%levelMeshes(currLevel)%numElements

   vtxLevel = currLevel - 1
   ! A vertex's parent is always a level above where the vertex forms another element
   ! Original mesh vertices have no parent, so should have 0 for both elementLevelId and level
   ! Level 1 vertices are on level 2
   ! For these original vertices
   ! create vertices
   do i = 1, size(vertices, dim=1)
      x = vertices(i,1)
      y = vertices(i,2)
      z = vertices(i,3)
      numNodes = numNodes + 1
      nodes(i) = createVertex(x, y, z, numNodes, domainId, vtxLevel)
   end do

   ! nodes(1) = createVertex(0.0d1, 0.0d1, 2.0d0, domainId, vtxLevel)
   ! numNodes = numNodes + 1
   ! nodes(2) = createVertex(0.1d1, -1.73d1, 0.0d1, numNodes, domainId, vtxLevel)
   ! ! nodes(2) = createVertex(0.0d1, 0.0d1, 0.0d1, domainId, vtxLevel)

   ! numNodes = numNodes + 1
   ! nodes(3) = createVertex(2.0d0, 0.0d1, 0.0d1, numNodes, domainId, vtxLevel)
   ! ! nodes(3) = createVertex(2.0d0, 0.0d1, 0.0d1, domainId, vtxLevel)

   ! numNodes = numNodes + 1
   ! nodes(4) = createVertex(0.1d1, -0.58d0, 0.163d1, numNodes, domainId, vtxLevel)
   ! ! Node 5 for 2nd tetrahedron
   ! numNodes = numNodes + 1
   ! nodes(5) = createVertex(0.267d1, -0.154d1, 1.09d0, numNodes, domainId, vtxLevel)
   ! nodes(4) = createVertex(0.0d1, 2.0d0, 0.0d1, domainId, vtxLevel)
   write(*,*) "Creating 1st element..."
   ! create first tetrahedron with connectivity -> [1, 2, 3, 4]
   elementNodes = (/nodes(1), nodes(2), nodes(3), nodes(4)/)
   numElements = numElements + 1
   numLevelElements = numLevelElements + 1

   ! For initial mesh, the parent level ID is 0
   tet(1) = createTetrahedron(elementNodes, numLevelElements, numElements, currLevel, 0)

   write(*,*) "Creating node for 2nd element..."

   ! nodes(5) = createVertex(0.0d1, 0.0d1, 2.0d0, domainId, vtxLevel)

   write(*,*) "Creating 2nd element..."
   ! create 2nd tet with connectivity -> [4, 2, 3, 5]
   elementNodes = (/nodes(4), nodes(2), nodes(3), nodes(5)/)
   numElements = numElements + 1
   numLevelElements = numLevelElements + 1
   tet(2) = createTetrahedron(elementNodes, numLevelElements, numElements, currLevel, 0)

   write(*,*) "Data of created Tetrahedrons"
   do i = 1, 2
      write(*,*) "For Tetrahedron ", i
      write(*,*) "Vertices: ", tet(i)%vertices
      write(*,*) "Element Level ID: ", tet(i)%elementLevelId
      write(*,*) "Element Global ID: ", tet(i)%globalId
      write(*,*) "Level: ", tet(i)%level
   end do

   write(*,*) "Adding tets as array to base_mesh..."
   ! allocate cell space in BaseMesh for tetrahedrons
   ! No need to allocate cell space when adding an array of tets only needed when adding
   ! one at a time
   call base_mesh%addTetrahedrons(currLevel, tet)

   write(*,*) "Tetrahedrons added successfully"
   write(*,*) "Base Mesh Data: "
   write(*,*) "Number of Elements: ", base_mesh%numElements
   write(*,*) "Number of Nodes: ", base_mesh%numVertices
   write(*,*) "For Nodes (V2pe): "
   do i = 1, base_mesh%numVertices
      write(*,*) "for node with global ID: ", i
      write(*,*) "Global ID on Vertex obj: ", base_mesh%v2pe(i)%globalId
      write(*,*) "Parent Element Level ID: ", base_mesh%v2pe(i)%elementLocalId
      write(*,*) "Parent Element Level: ", base_mesh%v2pe(i)%elementLevel
      write(*,*) "Vertex Local ID: ", base_mesh%v2pe(i)%localId
   end do
   write(*,*) "Number of Level Meshes: ", base_mesh%numLevels

   write(*,*) "------------------------------------------------"
   write(*,*) "Level 1 Mesh Data: "
   write(*,*) "Elements: "
   do i = 1, size(base_mesh%levelMeshes(currLevel)%elements)
      write(*,*) "For element Level ID ", i
      write(*,*) "Element global ID", base_mesh%levelMeshes(currLevel)%elements(i)
      write(*,*) "Element On base_mesh: "
      write(*,*) "confirming element level ID: ", &
         base_mesh%tetrahedrons(base_mesh%levelMeshes(currLevel)%elements(i))%elementLevelId
      write(*,*) "With edges: "
      do j = 1, 6
         write(*,*) "Edge: ", j
         write(*,*) base_mesh%tetrahedrons(base_mesh%levelMeshes(currLevel)%elements(i))%edges(j)
      end do
   end do
   write(*,*) "Number of elements: ", base_mesh%levelMeshes(currLevel)%numElements
   ! element connectivity no longer used
   ! write(*,*) "Element Connectivity: ", base_mesh%levelMeshes(currLevel)%elementConn

   write(*,*) "------------------------------------------------"

   write(*,*) "Building mesh connectivity for level 1..."
   ! build Mesh connectivity for this level 1 mesh
   call base_mesh%buildLevelMeshConnectivity(currLevel)

   write(*,*) "To check sibling half-facets (sibhfcs) array for level 1 mesh: "
!    write(*,*) "Level 1 Mesh sibhfcs: ", base_mesh%levelMeshes(currLevel)%sibhfcs

   do i = 1, size(base_mesh%levelMeshes(currLevel)%sibhfcs, dim=1)
      do j = 1, 4
         write(*,*) "Sibling Tuple for level element: ", i
         write(*,*) "for facet: ", j
         write(*,*) base_mesh%levelMeshes(currLevel)%sibhfcs(i,j)
      end do
   end do

   write(*,*) "**********************************************************"
   write(*,*) "To check vertex 2 half-facet(v2hf) array for level 1 mesh: "
   do i = 1, size(base_mesh%levelMeshes(currLevel)%v2hf)
      write(*,*) "for vertex ", i
      write(*,*) "Element ID: ", base_mesh%levelMeshes(currLevel)%v2hf(i)%elementLevelId
      write(*,*) "Facet ID: ", base_mesh%levelMeshes(currLevel)%v2hf(i)%localId
   end do
   write(*,*) "**********************************************************"

   ! To test creating a MeshEdge
   write(*,*) "Creating mesh edge..."
   edge = createMeshEdgeSorted(nodes(1)%globalID, nodes(2)%globalId)

   write(*,*) "Mesh Edge successfully created with vertices: "
   write(*,*) edge%vertices
   write(*,*) "Is edge Null: ", edge%isNullEdge()
   write(*,*) "Is edge equal: ", edge%isEdgeEqual(edge)

   ! write(*,*) "**********************************************************"
   ! write(*,*) "Testing FTL HashSet Int..."
   ! write(*,*) "Initializing FTL HashSet Int..."
   ! call intSet%New(5)
   ! call intSet%Insert(20)
   ! call intSet%Insert(11)
   ! call intSet%Insert(11)

   ! write(*,*) "Is intSet Empty: ", intSet%Empty()
   ! write(*,*) "Is 11 in intSet: ", intSet%Has(11)
   ! write(*,*) "Size of intSet: ", intSet%Size()
   ! write(*,*) "Elements in intSet: "
   ! intIter = intSet%Begin()
   ! do i = 1, intSet%Size()
   !    if ( intIter == intSet%End() ) exit
   !    write(*,*) intIter%value
   !    call intIter%Inc()
   !    ! write(*,*) intSet%End()%value
   ! end do

   ! call intSet2%New(5)
   ! call intSet2%Insert(30)
   ! call intSet2%Insert(11)
   ! call intSet2%Insert(10)
   ! call intSet2%Insert(1)

   ! write(*,*) "Is resultSet empty: ", resultSet%Empty()
   ! call resultSet%New(5)

   ! intIter = intSet%Begin()
   ! do i = 1, intSet%Size()
   !    if (intSet2%Has(intIter%value)) then
   !       call resultSet%Insert(intIter%value)
   !    end if
   !    call intIter%Inc()
   !    ! write(*,*) intSet%End()%value
   ! end do
   ! write(*,*) "To view result set..."
   ! ! resIter = resultSet%Begin()
   ! intIter = resultSet%Begin()

   ! do i = 1, resultSet%Size()
   !    !    if ( resIter == resultSet%End() ) exit
   !    !    write(*,*) resIter%value
   !    !    call resIter%Inc()
   !    if ( intIter == resultSet%End() ) exit
   !    write(*,*) intIter%value
   !    call intIter%Inc()
   ! end do
   write(*,*) "***********************************************************"
   write(*,*) "To check v2es of level 1 mesh: "
   do i = 1, size(base_mesh%levelMeshes(currLevel)%v2es)
      if ( base_mesh%levelMeshes(currLevel)%v2es(i)%vertexGlobalId /= -1 ) then
         write(*,*) "Element level IDs for Node: ", i
         write(*,*) base_mesh%levelMeshes(currLevel)%v2es(i)%v2elems
         write(*,*) "To confirm global node ID: ", base_mesh%levelMeshes(currLevel)%v2es(i)%vertexGlobalId
      end if
   end do
   write(*,*) "***********************************************************"
   write(*,*) "To check facet-attached neighbours of each element..."
   do i = 1, base_mesh%levelMeshes(currLevel)%numElements
      neighbours = base_mesh%getLevelMeshElementFacetAttachedNeighbours(i, currLevel)
      write(*,*) "Neighbours for element with levelId: ", i
      do j = 1, 4
         write(*,*) "For facet: ", j
         write(*,*) "Is there No Neighbour: ", neighbours(j)%isNullElement()
         if ( .not. neighbours(j)%isNullElement() ) then
            write(*,*) "The neighbour cell globalId: ", neighbours(j)%globalId
         end if
      end do
      write(*,*) "----------------------------------------------------------"
   end do

   ! allocate(intArr(6))
   ! write(*,*) "print intArr: ", intArr

   write(*,*) "**************************************************************"
   write(*,*) "Test binary search: "
   intArr = (/1, 3, 4, 7, 9, 20/)
   write(*,*) "index of 4: ", binarySearch(intArr, 4)
   write(*,*) "index of 8: ", binarySearch(intArr, 8)
   write(*,*) "index of 20: ", binarySearch(intArr, 20)
   write(*,*) "index of 1: ", binarySearch(intArr, 1)
   write(*,*) "**************************************************************"

   nodes(1) = base_mesh%getEdgeMidVertex(edge, base_mesh%tetrahedrons(1))
   write(*,*) "Test if element with level ID 1 has edge 1,2 and the midVertex: "
   write(*,*) "Is MidVertex null: ", nodes(1)%isNullVertex()
   write(*,*) "Does element have edge? ", base_mesh%tetrahedrons(1)%elementHasEdge(edge)

   write(*,*) "**************************************************************"
   write(*,*) "To test hanging nodes for L1: "
   hangingNodes = base_mesh%getLevelMeshHangingNodes(currLevel)
   write(*,*) "Size of hanging nodes: ", size(hangingNodes)

   write(*,*) "**************************************************************"
   write(*,*) "Test Refine level element for first element: "
   call base_mesh%refineLevelElement(currLevel, tet(1)%elementLevelId)
   write(*,*) "**************************************************************"
   write(*,*) "Test Refine level element for first child of first element: "
   nextLevel = currLevel + 1
   call base_mesh%refineLevelElement(nextLevel, base_mesh%levelMeshes(currLevel)%e2ce(1))

   write(*,*) "**************************************************************"
   write(*,*) "To check sibling half-facets (sibhfcs) array for level 2 mesh: "
!    write(*,*) "Level 1 Mesh sibhfcs: ", base_mesh%levelMeshes(currLevel)%sibhfcs

   do i = 1, size(base_mesh%levelMeshes(nextLevel)%sibhfcs, dim=1)
      do j = 1, 4
         write(*,*) "Sibling Tuple for level element: ", i
         write(*,*) "for facet: ", j
         write(*,*) base_mesh%levelMeshes(nextLevel)%sibhfcs(i,j)
      end do
   end do

   write(*,*) "**********************************************************"

   ! generate element files for every level
   do i = 1, base_mesh%numLevels
      write(*,*) "Attempting to generate element file for level: ", i
      call generateLevelMeshElementsFile(base_mesh, i)
      write(*,*) "File successfully generated!"
   end do
   write(*,*) "**********************************************************"
   write(*,*) "Attempting to generate nodes file..."
   call generateNodesFile(base_mesh)
   write(*,*) "File successfully generated!"
   write(*,*) "**********************************************************"
   ! generate hanging node files for every level
   do i = 1, base_mesh%numLevels
      write(*,*) "Attempting to generate hanging node file for level: ", i
      call generateLevelMeshHangingNodesFile(base_mesh, i)
      write(*,*) "File successfully generated!"
   end do
   write(*,*) "**********************************************************"
   ! generate mesh view files for every level
   do i = 1, base_mesh%numLevels
      write(*,*) "Attempting to generate mesh views file for level: ", i
      call generateLevelMeshViewFile(base_mesh, i)
      write(*,*) "File successfully generated!"
   end do
   write(*,*) "**********************************************************"
   ! generate node interpolation file
   write(*,*) "Attempting to generate nodes interpolation file"
   call generateNodesInterpolationFile(base_mesh)
   write(*,*) "File successfully generated!"
   write(*,*) "**********************************************************"
end program main
