program regular_tet
   use MeshLib
   use IOLib
   implicit none
   type(BaseMesh) :: base_mesh
   type(Vertex) :: nodes(10), elementNodes(4)
   type(Tetrahedron), allocatable :: tets(:)
   integer :: numNodes, currLevel, numElements, numLevelElements, i, domainId, vtxLevel
   integer, allocatable :: markedTets(:)
   integer :: elementConn(1,4)
   double precision :: vertices(4,3)
   double precision :: x, y, z
   ! Read entire nodes from file into an array that typically is number of nodes * 3 coords in size
!    vertices(1,:) = (/0.0000000000000000,0.0000000000000000, 2.0000000000000000/)
!    vertices(2,:) = (/2.0000000000000000, 0.0000000000000000, 0.0000000000000000/)
!    vertices(3,:) = (/0.0000000000000000, 0.0000000000000000, 0.0000000000000000/)
!    vertices(4,:) = (/1.0000000000000000, 2.0000000000000000, 1.0000000000000000/)
!    vertices(5,:) = (/1.0000000000000000, 0.0000000000000000, 1.0000000000000000/)
!    vertices(6,:) = (/0.0000000000000000, 0.0000000000000000, 1.0000000000000000/)
!    vertices(7,:) = (/0.50000000000000000, 1.0000000000000000, 1.5000000000000000/)
!    vertices(8,:) = (/1.0000000000000000,0.0000000000000000, 0.0000000000000000/)
!    vertices(9,:) = (/1.5000000000000000, 1.0000000000000000,0.50000000000000000/)
!    vertices(10,:) = (/0.50000000000000000, 1.0000000000000000,0.50000000000000000/)
   vertices(1,:) = (/0.0, 0.0, 0.0/)
   vertices(2,:) = (/1.0, -1.73, 0.0/)
   vertices(3,:) = (/2.0, 0.0, 0.0/)
   vertices(4,:) = (/1.0, -0.58, 1.63/)

   ! Read entire element connectivity file into array. You can actually decide to create the tetrahedrons
   ! directly
!    elementConn(1,:) = (/1, 5, 7, 8/)
!    elementConn(2,:) = (/5, 2, 6, 9/)
!    elementConn(3,:) = (/7, 6, 3, 10/)
!    elementConn(4,:) = (/8, 9, 10, 4/)
!    elementConn(5,:) = (/5, 7, 8, 9/)
!    elementConn(6,:) = (/5, 7, 6, 9/)
!    elementConn(7,:) = (/7, 8, 9, 10/)
!    elementConn(8,:) = (/7, 6, 9, 10/)
   elementConn(1,:) = (/1, 2, 3, 4/)
   ! The entire domain being meshed is taken to be on Level 0
   domainId = 0

   ! Allocate array to store actual tetrahedron elements
   allocate(tets(1))

   ! Initialize a BaseMesh: to set up Hi-AHF
   write(*,*) "initializing BaseMesh..."
   call base_mesh%initialize()
   write(*,*) "initializing LevelMesh..."
   ! initialize a levelMesh for level 1 - for the original, initial coarse mesh.
   ! available level meshes auto increment with each initialization, you typically only have to do this
   ! once, as the refinement progresses, levelMeshes are automatically added. Perhaps this should
   ! be hidden in the initialization of the base mesh
   call base_mesh%initializeLevelMesh()

   write(*,*) "successfully initialized BaseMesh and LevelMesh..."

   ! Create actual node objects (Vertex) struct.
   !  Allocate 1D array the size of the number of nodes for this
   write(*,*) "Creating nodes for 1st element..."

   ! For tracking. These parameters are required to create nodes and tetrahedrons
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

   ! Create actual tetrahedrons for the first level
   do i = 1, size(elementConn, dim=1)
      ! To track the global Id for the elements
      numElements = numElements + 1
      ! access the array of actual node objects via global node numbers specified in element connectivity
      elementNodes(1) = nodes(elementConn(i,1))
      elementNodes(2) = nodes(elementConn(i,2))
      elementNodes(3) = nodes(elementConn(i,3))
      elementNodes(4) = nodes(elementConn(i,4))
      ! Signature for this function:
      ! createTetrahedron(vertices, elementlevelid, globalid, level, parentelementlevelid)
      ! parentElementLevelId is 0 for elements of original mesh
      tets(i) = createTetrahedron(elementNodes, i, numElements, currLevel, 0)
   end do

   ! Add all the tetrahedrons just created to the base mesh Hi-AHF
   call base_mesh%addTetrahedrons(currLevel, tets)

   ! Build mesh connectivity for this original mesh. This must not be skipped
   call base_mesh%buildLevelMeshConnectivity(currLevel)
   ! to refine, read either the global Ids of elements you want to refine from file and save in a 1D
   ! array, or read the elementLevelIds and the element levels from file
   ! in this case we want to refine elements 2, 3, and 4, using the global Ids

   ! plan is to refine the first element into 8 tets, then refine children
   ! 3, 4, 5 of the parent element
   call base_mesh%refineLevelElement(currLevel, tets(1)%elementLevelId)
   ! at this point, we have 9 tets, 1 on level 1, 8 on level 2
   ! to refine elements 3, 4, 5 on level 2, we set the elementLevelIds
   markedTets = (/3, 4, 5/)

   ! loop over marked Tets and refine
   do i = 1, size(markedTets)
      ! get the element from the base mesh if using globalIds
      !   markedTet = base_mesh%tetrahedrons(markedTets(i))
      ! get the level and elementLevelId
      !   call base_mesh%refineLevelElement(markedTet%level, markedTet%elementLevelId)
      call base_mesh%refineLevelElement(2, markedTets(i))
      ! You can easily generate files here as this progresses but they will keep being rewritten
      ! maybe an argument to them so they can generate different files?
   end do
   ! generate nodes file
   call generateNodesFile(base_mesh)

   ! for the per level files, generate by looping over the available levels
   do i = 1, base_mesh%numLevels
      ! generate elements file
      call generateLevelMeshElementsFile(base_mesh, i)
      ! generate hanging nodes file
      call generateLevelMeshHangingNodesFile(base_mesh, i)
      ! generate visualization file
      call generateLevelMeshViewFile(base_mesh, i)
   end do

   ! generate global nodes interpolation file
   call generateNodesInterpolationFile(base_mesh)
   ! compile helpers.f90, MeshLib.f90, and IOLib.f90 modules into your program, and run
   ! gfortran -g -fcheck=all ./MeshLib.f90 ./helpers.f90 ./IOLib.f90 ./regular_tet.f90 -o regular_tet
end program regular_tet
