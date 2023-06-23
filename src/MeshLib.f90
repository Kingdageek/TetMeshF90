module MeshLib
   ! use stdlib_sorting, only: ord_sort
   ! use ftlHashSetIntModule
   use helpers
   implicit none

   ! ONLY FOR IRREGULAR REFINEMENT
   ! Edge entity to help with refinement and to avoid duplicate nodes
   ! Will typically only be created on the fly
   type :: MeshEdge
      ! Global IDs of the 2 connecting vertices (preferably sorted in ASC)
      ! Default values for the vertices to ensure we can check for a NULL edge
      integer :: vertices(2) = (/0,0/)

   contains
      procedure :: isNullEdge, isEdgeEqual
   end type MeshEdge

   type :: Vertex
      ! Coordinates
      double precision :: x, y, z
      ! Global Node Index. NULL vertex will have globalId = -1
      integer :: globalId = -1

      ! A vertex can belong to multiple elements
      ! Id of parent element on concerned level that vertex belongs to
      !!! These ppties below are tighly coupled to an element and a vertex can be shared by several elements
      ! Id of one of its parent cells/tetrahedrons/elements
      ! Needed in vertex to parent mapping (v2pe) in BaseMesh. Level 1 vertices are actually on Level 2
      integer :: elementLocalId = -1
      ! intended to be the local ID of this vertex on the current element it bounds
      integer :: localId = -1
      ! local ID of vertex on the element it was interpolated from i.e. its parent, follows the CGNS
      ! It would also be the same as the local ID of the refined Edge i.e.
      ! the local ID of the edge in the parent element with `elementLocalId`
      ! This will now be set during refinement so there'd be less stress when creating
      ! original mesh vertices
      integer :: elementVertexId = -1
      ! level of parent this vertex was obtained from i.e. element it was
      ! interpolated from
      integer :: elementLevel = -1

      ! nodes from the parent element that this node was interpolated from
      ! this node is a midnode of the connecting nodes of this edge
      ! if this edge is NULL after processing => this NODE is an original domain vertex
      type(MeshEdge) :: parentEdge

   contains
      procedure :: isNullVertex
   end type Vertex

   ! half-face of tetrahedron (3D) - triangle
   type :: HalfFacet
      ! connecting nodes ordered following the
      ! CFD General Notation System (CGNS)
      type(Vertex) :: vertices(3)
      ! Id of element that owns this HalfFacet on the mesh of the level concerned
      integer :: elementLevelId = -1
      ! local Id of this facet on the element, follows CGNS too
      integer :: localId = -1

   contains
      procedure :: isHfEqual, isNullHf
   end type HalfFacet

   ! Actual 3D element
   type :: Tetrahedron
      ! Connectivity of Vertex
      type(Vertex) :: vertices(4)
      ! HalfFacets that form this element
      type(HalfFacet) :: facets(4)
      ! 6 edges for this tetrahedron, with sorted global node numbers
      type(MeshEdge) :: edges(6)
      ! current refinement level
      integer :: level = -1
      ! global element index
      integer :: globalId = -1
      ! ID of this element on its level
      integer :: elementLevelId = -1

      ! ID of this element's parent on level -1
      integer :: parentElementLevelId = -1

      ! Flag to check if element has been refined
      ! An element is said to be 'ACTIVE' if it hasn't been refined
      ! every created element should be active by default
      logical :: isActive = .true.
   contains
      procedure :: isNullElement, elementHasEdge
   end type Tetrahedron

   type :: MeshConstants
      real :: v2peIncrement = 0.2
   end type MeshConstants

   type :: SiblingHFTuple
      ! Each cell is a tuple of form:
      ! <sibling half facet element level Id, sibhfc local Id on element>
      ! default values of -1 if this tuple is NOT SET in sibhfcs
      ! values set to 0 if facet is on boundary, i.e. does not have a sibling facet
      ! integer :: sibhfElemLevelId = -1
      ! integer :: sibhfLocalId = -1

      ! Let every sibhfTuple represent that half-facet is all boundary by default
      ! we change this in sibhfcs only for facets that have actual siblings
      integer :: sibhfElemLevelId = 0
      integer :: sibhfLocalId = 0

   contains
      procedure :: isNullSibHFTuple
   end type SiblingHFTuple

   ! For the intermediate structure v2hfs (mapping of vertex to its incident
   ! half-facets in which the vertex has the largest ID)
   type :: Vertex2HalfFacets
      integer :: vertexGlobalId = -1
      ! Array of every HF incident on vertex with `vertexGlobalId`
      type(HalfFacet), allocatable :: hfs(:)
   end type Vertex2HalfFacets

   ! for the intermediate structure v2adj (mapping of vertex to its adjacent vertices
   ! in each of the incident half-facets in v2hfs)
   type :: VertexHalfFacet
      integer :: vertexGlobalId = -1
      integer :: elementGlobalId = -1
      integer :: elementLevelId = -1
      ! local Id of this facet in element with `elementGlobalId`
      integer :: facetId = -1
   end type VertexHalfFacet

   ! in v2adj, vertex => [[adjacent vertices on facet1], [adjacent vertices on facet2]]
   ! [
   !  v0 => [[VHF1_1, VHF1_2, VHF1_3], [VHF2_1, VHF2_2, VHF2_3]],
   !  v1 => [[V1HF1_1, V1HF1_2, V1HF1_3]]
   !]
   type :: Vertex2VertexHalfFacets
      integer :: vertexGlobalId = -1
      type(VertexHalfFacet), allocatable :: vhf(:,:)
   end type Vertex2VertexHalfFacets

   type :: Vertex2Elements
      integer :: vertexGlobalId = -1
      integer, allocatable :: v2elems(:)
   end type Vertex2Elements
   ! Mesh representation for a level using Array-based Half-Facets
   type :: LevelMesh
      ! connectivity matrix for each element of this level. Just topological info
      ! shape: numOfELementsInLevel x 4 vertice global Ids per element: ((1, 2, 3, 4), (11, 12, 13, 14))
      ! we can just save the globalIds of the elements here and literally fetch them from: tet= base_mesh%tetrahedrons(globalId)
      ! For connectivity: tet%vertices(1)%globalId
      integer, allocatable :: elementConn(:,:)
      ! to store the global IDs of all elements on this level for retrieval from the
      ! base_mesh
      integer, allocatable :: elements(:)
      ! sibling half-facets: AHF data for adjacency queries
      ! index is levelIndex for elements,
      ! shape: number of elements on level x 4 facets per element
      ! Each cell is a tuple of form: <sibling half facet element level Id, sibhfc local Id on element>
      ! Rough example: [[0, <2,3>, 0, 0], [<3,3>, <1,1>, 0, 0]]
      ! [[NULL, SHFT<shfElemLevelId:2, shfLocalId3>, NULL, NULL]]
      type(SiblingHFTuple), allocatable :: sibhfcs(:,:)

      ! vertex to half-facets, anchor to locate vertex
      ! to help determine border vertices easily
      ! vertex to an incident half-facet
      type(HalfFacet), allocatable :: v2hf(:)

      ! vertex to elements array. a map of v -> hashset of element global IDs
      ! It's a map of each vertex to every element on a level incident on it
      ! To help find elements on a level attached to a specified edge
      ! No hashset support anymore, we use an integer array: 2D array now
      ! type(ftlHashSetInt), allocatable :: v2es(:)
      type(Vertex2Elements), allocatable :: v2es(:)
      ! Element to Parent Element array for new elements on this level
      ! e2pe(elementLevelId) -> parentElementLevelId
      integer, allocatable :: e2pe(:)
      ! element to child element mapping. Maps to index of first child on its level
      ! other children are contiguous to this - array
      ! for each element onlevel, point to their 1st child on next level
      integer, allocatable :: e2ce(:)

      ! current number of elements on this level
      integer :: numElements

      ! level mesh belongs to
      integer :: level

   end type LevelMesh

   ! type :: VertexToParent
   !    integer ::
   ! end type VertexToParent

   ! Hierarchical Array-based Half-Facet(HAHF)
   type :: BaseMesh
      ! Statistics
      integer :: numVertices, numElements, numLevels
      ! mesh hierarchy. All meshes ordered with increasing level
      type(LevelMesh), allocatable :: levelMeshes(:)
      ! To track all elements and to be able to find element for refinement in O(1) time
      type(Tetrahedron), allocatable :: tetrahedrons(:)
      ! Nodes ordered according to global node index. Vertex to Parent Element
      ! stores tuples of form: <level, elementId, localId>. Stores new vertices.
      ! Original vertices are exempt i.e. Level 1 vertices are new vertices on L1 obtained by interpolation -> midpoint vertices
      type(Vertex), allocatable :: v2pe(:)

      ! All vertices including original vertices. can be used for the v2pe mapping too
      ! The index is the node global index mapping to <level, elementId, localId>, for the original vertices <0, 0, localId>
      ! element 0 is the full domain, level 0 is level for full domain
      type(Vertex), allocatable :: vertices(:)

   contains
      procedure :: addTetrahedron, addTetrahedrons, initialize, initializeLevelMesh
      procedure :: allocateCellSpace, mapElementVerticesToParent, expandV2pe, buildLevelMeshConnectivity
      procedure :: buildLevelMeshV2hf, buildLevelMeshSibhfcs, buildLevelMeshV2es, getLevelMeshElementFacetAttachedNeighbours
      procedure :: getLevelElementsAttachedToEdge, getLevelElementIdsAttachedToVertex, getEdgeMidVertex
      procedure :: buildLevelMeshInterLevelConnectivity, getLevelMeshHangingNodes
      procedure :: refineLevelElement, createMidVertexOnEdge, localRegularRefinement
      procedure :: getLevelMeshNumNodes, getLevelMeshElementNeighbours, getElementLargeNeighbours
   end type BaseMesh

contains

   !!!!!! FUNCTIONS !!!!!!
   function createVertex(x, y, z, globalId, parentLevelId, parentLevel) result(v)
      implicit none
      double precision, intent(in) :: x, y, z
      integer, intent(in) :: globalId, parentLevelId, parentLevel
      type(Vertex) :: v
      v%x = x
      v%y = y
      v%z = z
      v%globalId = globalId
      ! this has to be set by the parent that this vertex bounds,
      ! but then this will change with each parent that shares this vertex that's created
      ! OR anytime a gro is created, before the element is created, check:
      ! levelMesh%numElements and add 1 to know the parentElementLevelId and pass that
      v%elementLocalId = parentLevelId
      v%elementLevel = parentLevel
      ! local cannonical Id of this vertex on the Parent it was interpolated from, also
      ! the ID of the edge this vertex is on. Let this be set during refinement
      ! v%elementVertexId = parentVertexId
      ! let local ID be set when v2pe is set
   end function createVertex

   ! creates an edge from the vertice global ID arguments
   ! sorts the global IDs in ASC order in edge%vertices
   function createMeshEdgeSorted(v1, v2) result(edge)
      implicit none
      integer, intent(in) :: v1, v2
      type(MeshEdge) :: edge
      integer :: tempMax, vertices(2)
      vertices = (/v1, v2/)
      tempMax = v1
      if ( tempMax > v2 ) then
         ! swap
         vertices(1) = v2
         vertices(2) = tempMax
      end if
      edge%vertices = vertices
   end function createMeshEdgeSorted

   function createHalfFacet(vertices, elementLevelId, localId) result(hf)
      implicit none
      type(Vertex), intent(in) :: vertices(3)
      integer, intent(in) :: elementLevelId, localId
      type(HalfFacet) :: hf

      hf%vertices = vertices
      hf%elementLevelId = elementLevelId
      hf%localId = localId
   end function createHalfFacet

   function createTetrahedron(vertices, elementLevelId, globalId, level, parentElementLevelId) result(tet)
      implicit none
      type(Vertex), intent(in) :: vertices(4)
      integer, intent(in) :: elementLevelId, globalId, level, parentElementLevelId
      type(Tetrahedron) :: tet
      type(HalfFacet) :: facets(4)
      type(MeshEdge) :: edges(6)

      ! vertices must be ordered following CGNS
      tet%vertices = vertices
      tet%elementLevelId = elementLevelId
      tet%globalId = globalId
      tet%level = level
      tet%parentElementLevelId = parentElementLevelId
      ! create half-facets for tet order by CGNS
      ! counter-clockwise to outward normal from the face, starting from left angle
      facets(1) = createHalfFacet((/vertices(1), vertices(3), vertices(2)/), elementLevelId, 1)
      facets(2) = createHalfFacet((/vertices(1), vertices(2), vertices(4)/), elementLevelId, 2)
      facets(3) = createHalfFacet((/vertices(2), vertices(3), vertices(4)/), elementLevelId, 3)
      facets(4) = createHalfFacet((/vertices(3), vertices(1), vertices(4)/), elementLevelId, 4)

      tet%facets = facets

      ! create edges for tet order by CGNS
      edges(1) = createMeshEdgeSorted(vertices(1)%globalId, vertices(2)%globalId)
      edges(2) = createMeshEdgeSorted(vertices(2)%globalId, vertices(3)%globalId)
      edges(3) = createMeshEdgeSorted(vertices(3)%globalId, vertices(1)%globalId)
      edges(4) = createMeshEdgeSorted(vertices(1)%globalId, vertices(4)%globalId)
      edges(5) = createMeshEdgeSorted(vertices(2)%globalId, vertices(4)%globalId)
      edges(6) = createMeshEdgeSorted(vertices(3)%globalId, vertices(4)%globalId)
      tet%edges = edges
   end function createTetrahedron

   function isNullVertex(vtx) result(isNull)
      implicit none
      class(Vertex), intent(in) :: vtx
      logical :: isNull
      isNull = .false.
      if ( vtx%globalId == -1 ) isNull = .true.
   end function isNullVertex

   function isNullEdge(edge) result(isNull)
      implicit none
      class(MeshEdge), intent(in) :: edge
      logical :: isNull
      integer :: v1, v2
      isNull = .false.
      v1 = edge%vertices(1)
      v2 = edge%vertices(2)

      if ( v1 == 0 .and. v2 == 0 ) then
         isNull = .true.
      end if
   end function isNullEdge

   ! checks if two edges are the same, i.e. are connected by the same 2 vertices
   ! vertices MUST be sorted in ASC order
   function isEdgeEqual(edge, edge2) result(isEqual)
      implicit none
      class(MeshEdge), intent(in) :: edge
      type(MeshEdge), intent(in) :: edge2
      logical :: isEqual
      isEqual = .false.

      if ( edge%vertices(1) == edge2%vertices(1) .and. edge%vertices(2) == edge2%vertices(2)) then
         isEqual = .true.
      end if
   end function isEdgeEqual

   function isNullSibHFTuple(sibHft) result(isNull)
      class(SiblingHFTuple), intent(in) :: sibHft
      logical :: isNull
      isNull = .false.
      if ( sibHft%sibhfElemLevelId == 0 .or. sibHft%sibhfLocalId == 0 ) then
         isNull = .true.
      end if
   end function isNullSibHFTuple

   function getHfVertexWithLargestId(hf) result(vtxMax)
      implicit none
      type(HalfFacet), intent(in) :: hf
      type(Vertex) :: vtxMax
      integer :: i

      ! set max global ID to -1, to make max a null vertex
      vtxMax%globalId = -1

      ! over the vertices of this facet. Finding the MAX
      do i = 1, 3
         if ( vtxMax%globalId < hf%vertices(i)%globalId ) then
            vtxMax = hf%vertices(i)
         end if
      end do
   end function getHfVertexWithLargestId

   function isNullHf(hf) result(isNull)
      implicit none
      class(HalfFacet), intent(in) :: hf
      logical :: isNull
      isNull = .false.
      if ( hf%elementLevelId == -1 .or. hf%localId == -1 ) then
         isNull = .true.
      end if
   end function isNullHf

   function isSiblingHFTupleSet(shfTuple) result(isSet)
      implicit none
      type(SiblingHFTuple), intent(in) :: shfTuple
      logical :: isSet
      isSet = .false.
      if (shfTuple%sibhfElemLevelId > 0 .and. shfTuple%sibhfLocalId > 0) then
         isSet = .true.
      end if
   end function isSiblingHFTupleSet

   function isSiblingHFTupleOnBoundary(shfTuple) result(isOnBoundary)
      implicit none
      type(SiblingHFTuple), intent(in) :: shfTuple
      logical :: isOnBoundary
      isOnBoundary = .false.
      if (shfTuple%sibhfElemLevelId == 0 .and. shfTuple%sibhfLocalId == 0) then
         isOnBoundary = .true.
      end if
   end function isSiblingHFTupleOnBoundary

   function getAdjacentVerticesToMaxVtxInFacet(v2adj, vtxMaxGlobalId, hf) result (vhfs)
      implicit none
      type(Vertex2VertexHalfFacets), allocatable, intent(in) :: v2adj(:)
      integer, intent(in) :: vtxMaxGlobalId
      type(HalfFacet), intent(in) :: hf
      ! always 2 adjacent vertices to max vertex per face
      type(VertexHalfFacet) :: vhfs(2)
      ! 2D: number of half-facets adjacent to vtxMax x 2 possible adj. vertices per facet
      TYPE(VertexHalfFacet), ALLOCATABLE, DIMENSION(:,:) :: vhf
      integer :: i

      ! To check for adjacent vertices on `hf` where vtx with `vtMaxGlobalId` is max
      vhf = v2adj(vtxMaxGlobalId)%vhf
      do i = 1, size(vhf)
         ! check first cell to find facet `hf`, if there's a match, the entire row
         ! is a set of vertices adj. to vtxMax on facet `hf`
         if ( vhf(i,1)%elementLevelId == hf%elementLevelId .and. &
            vhf(i,1)%facetId == hf%localId ) then
            vhfs = vhf(i,:)
            ! adjacent vertices found, exit loop
            exit
         end if
      end do
   end function getAdjacentVerticesToMaxVtxInFacet

   function isHfEqual(hf1, hf2) result(isEqual)
      implicit none
      class(HalfFacet), intent(in) :: hf1
      type(HalfFacet), intent(in) :: hf2
      logical :: isEqual
      isEqual = .false.
      if ( (hf1%localId == hf2%localId) .and. &
         (hf1%elementLevelId == hf2%elementLevelId) ) then
         isEqual = .true.
      end if
   end function isHfEqual

   function isSameAdjacentVertices(adjVhfs, tempAdjVhfs) result (isSame)
      implicit none
      TYPE(VertexHalfFacet), DIMENSION(2), intent(in) :: adjVhfs
      TYPE(VertexHalfFacet), DIMENSION(2), intent(in) :: tempAdjVhfs
      logical :: isSame
      ! isSame = .false.
      ! to be true, facetId and elementGlobalId at 1st position of both adjVhfs & tempAdjVhfs
      ! must match, similar for the 2nd position. To be sure the adjacent vertices are the
      ! same OR we can just compare the globalId of the Vertices
      ! isSame = (adjVhfs(1)%facetId == tempAdjVhfs(1)%facetId .and. &
      !    adjVhfs(1)%elementGlobalId == tempAdjVhfs(1)%elementGlobalId) .and. &
      !    (adjVhfs(1)%facetId == tempAdjVhfs(1)%facetId .and. &
      !    adjVhfs(1)%elementGlobalId == tempAdjVhfs(1)%elementGlobalId)
      isSame = (adjVhfs(1)%vertexGlobalId == tempAdjVhfs(1)%vertexGlobalId) .and. &
         (adjVhfs(2)%vertexGlobalId == tempAdjVhfs(2)%vertexGlobalId)
   end function isSameAdjacentVertices

   function getHalfFacetsWithSameAdjacentVertices(v2hfs, vtxMaxGlobalId, v2adj, adjVhfs) &
      result(hfs)
      implicit none
      TYPE(Vertex2HalfFacets), ALLOCATABLE, DIMENSION(:), intent(in) :: v2hfs
      integer, intent(in) :: vtxMaxGlobalId
      TYPE(Vertex2VertexHalfFacets), ALLOCATABLE, DIMENSION(:), intent(in) :: v2adj
      ! set of adjacent vertices (mapped to the facets) to compare against
      TYPE(VertexHalfFacet), DIMENSION(2), intent(in) :: adjVhfs
      ! array of matching half-facets to be returned
      type(HalfFacet), allocatable :: hfs(:)
      type(HalfFacet), allocatable :: tempHfs(:)
      integer :: i, j, k
      TYPE(VertexHalfFacet), DIMENSION(2) :: tempAdjVhfs
      type(HalfFacet) :: hf

      allocate(hfs(0))
      ! to track added half-facets
      k = 1
      ! for each half-facet, hf, in v2hfs(vtxMaxGlobalId)%hfs
      ! get v2adj(vtxMax, hf)
      ! compare to adjVhfs, if similar, add to hfs array
      do i = 1, size(v2hfs(vtxMaxGlobalId)%hfs)
         hf = v2hfs(vtxMaxGlobalId)%hfs(i)
         tempAdjVhfs = getAdjacentVerticesToMaxVtxInFacet(v2adj, vtxMaxGlobalId, hf)
         ! they should be sorted, so we can compare positional indices directly
         if ( isSameAdjacentVertices(adjVhfs, tempAdjVhfs) ) then
            ! expand array to contain next element
            if (allocated(tempHfs)) deallocate(tempHfs)
            allocate(tempHfs(size(hfs) + 1))
            ! copy to hfs
            do j = 1, size(hfs)
               tempHfs(j) = hfs(j)
            end do
            tempHfs(k) = hf
            deallocate(hfs)
            hfs = tempHfs
            k = k + 1
         end if
      end do

   end function getHalfFacetsWithSameAdjacentVertices

   function getLevelMeshElementFacetAttachedNeighbours(base_mesh, elementLevelId, meshLevel) &
      result(neighbours)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: elementLevelId, meshLevel
      type(LevelMesh) :: mesh
      ! only 4 facets per element -> max of 4 neighbours attached to a facet
      ! Null tetrahedrons may be included if the facet checked is a border facet
      type(Tetrahedron) :: neighbours(4)
      integer :: i
      type(SiblingHFTuple) :: sibhfTuple

      mesh = base_mesh%levelMeshes(meshLevel)

      ! over the elements facets
      do i = 1, 4
         sibhfTuple = mesh%sibhfcs(elementLevelId, i)
         if ( .not. sibhfTuple%isNullSibHFTuple() ) then
            ! mesh%elements(elemId) returns the globalId of the element in base_mesh
            neighbours(i) = base_mesh%tetrahedrons(mesh%elements(sibhfTuple%sibhfElemLevelId))
         end if
      end do
   end function getLevelMeshElementFacetAttachedNeighbours

   function getLevelElementIdsAttachedToVertex(base_mesh, meshLevel, vtxGlobalId) result(elements)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel, vtxGlobalId
      integer, allocatable :: elements(:)
      type(LevelMesh) :: mesh

      mesh = base_mesh%levelMeshes(meshLevel)
      elements = mesh%v2es(vtxGlobalId)%v2elems
   end function getLevelElementIdsAttachedToVertex

   function getLevelElementsAttachedToEdge(base_mesh, meshLevel, edge) result(elements)
      class(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel
      type(MeshEdge), intent(in) :: edge
      type(Tetrahedron), allocatable :: elements(:), tempElems(:)
      type(LevelMesh) :: mesh
      integer, allocatable :: v1Cells(:), v2Cells(:)
      integer :: i, j, elemLevelId, k
      type(Tetrahedron) :: tet

      mesh = base_mesh%levelMeshes(meshLevel)
      ! get level elements attached to edge%v1 and the ones attached to edge%v2
      ! find the intersection of these elements
      ! should have at least one element, the element for which this check happens
      ! calling for edges or vertices which don't form an element should have weird effects
      v1Cells = base_mesh%getLevelElementIdsAttachedToVertex(meshLevel, edge%vertices(1))
      v2Cells = base_mesh%getLevelElementIdsAttachedToVertex(meshLevel, edge%vertices(2))

      allocate(elements(0))
      ! find intersection
      ! a set will allow for constant time search. No set here but we take advantage of the
      ! idea that the cells will be in increasing order cos mesh%v2es starts comparison
      ! from 1 to maxNumberElements, => we can use a binary search to do checks in O(log(n)) time
      do i = 1, size(v1Cells)
         elemLevelId = v1Cells(i)
         ! search v2Cells
         ! j -> position of elemLevelId in v2Cells
         j = binarySearch(v2Cells, elemLevelId)
         if ( j /= -1 ) then
            ! write(*,*) "for element level ID: ", elemLevelId
            ! write(*,*) "V2Cells match is : ", v2Cells(j)
            tet = base_mesh%tetrahedrons(mesh%elements(elemLevelId))
            ! to add to the elem
            if (allocated(tempElems)) deallocate(tempElems)
            allocate(tempElems(size(elements) + 1))
            ! copy to tempElems
            do k = 1, size(elements)
               tempElems(k) = elements(k)
            end do
            ! add new element last
            tempElems(size(elements) + 1) = tet
            if (allocated(elements)) deallocate(elements)
            elements = tempElems
         end if
      end do
   end function getLevelElementsAttachedToEdge

   ! to get both edge and facet-attached neighbours
   function getLevelMeshElementNeighbours(base_mesh, elementLevelId, meshLevel) result(neighbours)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: elementLevelId, meshLevel
      type(Tetrahedron), allocatable :: tempNeighbours(:), neighbours(:), edgeNeighbours(:)
      ! type(Tetrahedron) :: facetNeighbours(4)
      type(LevelMesh) :: mesh
      integer :: i, j, levelId, numNeighbours
      type(Tetrahedron) :: checkedTet
      mesh = base_mesh%levelMeshes(meshLevel)

      ! allocate tempNeighbours to the numElements of this level Mesh - using an array as a map
      ! Necessary to remove duplicate elements quickly, to avoid comparing returned elements cos
      ! diff. edges and facets can have similar elements incidented on them
      allocate(tempNeighbours(mesh%numElements))
      ! There's actually no need to get the facet Neighbours cos Edge Neighbours will ALWAYS
      ! have every facet Neighbour Included
      ! for facet neighbours
      ! facetNeighbours = base_mesh%getLevelMeshElementFacetAttachedNeighbours(elementLevelId, meshLevel)
      ! ! facetNeighbours include NULL elements for boundary facets
      ! do i = 1, 4
      !    if ( .not. facetNeighbours(i)%isNullElement() ) then
      !       levelId = facetNeighbours(i)%elementLevelId
      !       tempNeighbours(levelId) = facetNeighbours(i)
      !    end if
      ! end do

      checkedTet = base_mesh%tetrahedrons(mesh%elements(elementLevelId))
      ! for edge neighbours - 6 edges per tet
      do i = 1, 6
         ! get every cell incidented on the edge, excluding the current cell
         edgeNeighbours = base_mesh%getLevelElementsAttachedToEdge(meshLevel, checkedTet%edges(i))
         ! Look through these and either simply update (if there was an element there already)
         ! tempNeighbours or if it's new slot it in
         do j = 1, size(edgeNeighbours)
            ! edgeNeighbours should not contain any NULL element.
            ! The checked TET is also included in cells incident on said edge so, check
            ! to avoid including it too
            levelId = edgeNeighbours(j)%elementLevelId
            if ( levelId /= checkedTet%elementLevelId ) then
               tempNeighbours(levelId) = edgeNeighbours(j)
            end if
         end do
      end do

      ! To avoid returning NULL elements, we count the actual neighbours found to initialize
      ! the returned array - neighbours
      numNeighbours = 0
      do i = 1, mesh%numElements
         if ( .not. tempNeighbours(i)%isNullElement() ) then
            numNeighbours = numNeighbours + 1
         end if
      end do

      allocate(neighbours(numNeighbours))
      ! to track position in neighbours
      j = 0
      do i = 1, mesh%numElements
         if ( .not. tempNeighbours(i)%isNullElement() ) then
            j = j+1
            neighbours(j) = tempNeighbours(i)
         end if
      end do
      write(*,*) "is j same as size of neighbours: ", j == size(neighbours)
   end function getLevelMeshElementNeighbours

   function getElementLargeNeighbours(base_mesh, tet) result(largeNeighbours)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      type(Tetrahedron), intent(in) :: tet
      type(Tetrahedron), allocatable :: largeNeighbours(:), tempNeighbours(:), edgeNeighbours(:)
      integer :: parentLevel, i, j, levelId, numNeighbours
      type(LevelMesh) :: mesh
      type(Vertex) :: vtx

      ! allocate largeNeighbours to size 0 to return that when it checks
      allocate(largeNeighbours(0))
      ! There'd be no largeNeighbours for elements on Level 1
      if (tet%level < 2) return
      parentLevel = tet%level - 1
      mesh = base_mesh%levelMeshes(parentLevel)
      ! allocate tempNeighbours to numElements on parentLevel mesh to help with avoiding duplicates
      ! and null values
      allocate(tempNeighbours(mesh%numElements))
      ! get parent edge for every vertex in tet
      ! over the vertices in tet
      do i = 1, 4
         vtx = tet%vertices(i)
         ! We need to check only Vertices INTERPOLATED from the parentLevel
         ! cos an element can have vertices created higher up the hierarchy too
         ! but there must be vertices directly interpolated from its parentElement
         if ( vtx%elementLevel /= parentLevel) cycle
         ! get elements incident to those edges
         ! get neighbours for the parentEdge of this vertex in the parent mesh
         ! write(*,*) "ParentLevel: ", parentLevel
         ! write(*,*) "ParentEdge vertices: ", vtx%parentEdge%vertices
         edgeNeighbours = base_mesh%getLevelElementsAttachedToEdge(parentLevel, vtx%parentEdge)
         do j = 1, size(edgeNeighbours)
            ! edgeNeighbours should not contain any NULL element.
            ! The parent of this TET is also included in cells incident on said edge so, check
            ! to avoid including it too.
            ! We can't count numNeighbours here cos diff. edges can have the same elements incident on them
            ! i.e. tempNeighbours(element) can be updated and re-updated, but there'd be just one eventually
            levelId = edgeNeighbours(j)%elementLevelId
            if ( levelId /= tet%parentElementLevelId ) then
               tempNeighbours(levelId) = edgeNeighbours(j)
            end if
         end do
      end do

      ! To avoid returning NULL elements, we count the actual neighbours found to initialize
      ! the returned array - neighbours
      numNeighbours = 0
      do i = 1, mesh%numElements
         if ( .not. tempNeighbours(i)%isNullElement() ) then
            numNeighbours = numNeighbours + 1
         end if
      end do
      if (allocated(largeNeighbours)) deallocate(largeNeighbours)
      allocate(largeNeighbours(numNeighbours))
      ! to track position in neighbours
      j = 0
      do i = 1, mesh%numElements
         if ( .not. tempNeighbours(i)%isNullElement() ) then
            j = j+1
            largeNeighbours(j) = tempNeighbours(i)
         end if
      end do
   end function getElementLargeNeighbours

   logical function elementHasEdge(tet, edge)
      implicit none
      class(Tetrahedron), intent(in) :: tet
      class(MeshEdge), intent(in) :: edge
      integer :: i
      elementHasEdge = .false.
      do i = 1, 6
         if ( tet%edges(i)%isEdgeEqual(edge) ) then
            elementHasEdge = .true.
            exit
         end if
      end do
   end function elementHasEdge

   function getEdgeMidVertex(base_mesh, edge, tet) result(midVertex)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      type(MeshEdge), intent(in) :: edge
      type(Tetrahedron), intent(in) :: tet
      type(Vertex) :: midVertex
      integer :: i, j, child1LevelId
      type(LevelMesh) :: mesh, nextMesh
      type(Tetrahedron) :: childTet

      ! if tet is a null element, return
      if ( tet%isNullElement() ) return
      mesh = base_mesh%levelMeshes(tet%level)
      ! check the children of tet for a node whose parent node is `edge`
      ! NULL vertex if tet is active => not refined or
      ! current level e2ce array is of size 0 or the next level doesn't exist
      write(*,*) "Check if tet is active: ", tet%isActive
      write(*,*) "Check if tet has edge: ", tet%elementHasEdge(edge)
      if ( .not. tet%isActive .and. tet%elementHasEdge(edge)) then
         ! tet's been refined, children should exist
         child1LevelId = mesh%e2ce(tet%elementLevelId)
         ! check the children for a node whose parent is edge
         nextMesh = base_mesh%levelMeshes(tet%level + 1)
         ! the eight children are stored in order
         do i = child1LevelId, (child1LevelId + 7)
            childTet = base_mesh%tetrahedrons(nextMesh%elements(i))
            ! over the vertices of child tet
            do j = 1, 4
               if ( childTet%vertices(j)%parentEdge%isEdgeEqual(edge) ) then
                  midVertex = childTet%vertices(j)
                  return
               end if
            end do
         end do
      end if

   end function getEdgeMidVertex

   logical function isNullElement(tet)
      implicit none
      class(Tetrahedron), intent(in) :: tet
      isNullElement = .false.
      if ( tet%elementLevelId == -1 .or. tet%level == -1 .or. tet%globalId == -1 ) then
         isNullElement = .true.
      end if
   end function isNullElement

   function getLevelMeshHangingNodes(base_mesh, meshLevel) result(hangingNodes)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel
      type(Vertex), allocatable :: hangingNodes(:), tempNodes(:)
      ! type(LevelMesh) :: mesh
      integer :: i,j,k, numHangingNodes, parentLevel
      type(Vertex) :: vtx
      type(Tetrahedron), allocatable :: tets(:)
      ! allocate it's size to 0, so it returns an empty array if no hanging node on level
      allocate(hangingNodes(0))
      numHangingNodes = 0
      parentLevel = (meshLevel - 1)
      ! No hanging Nodes on original mesh, we expect an initial conformal and manifold mesh
      if ( meshLevel > 1 ) then
         ! allocate size for all vertices, to minimize all the copying and size change
         ! You'd extract hanging nodes later to hanging nodes by checking if the vertex is null
         allocate(tempNodes(base_mesh%numVertices))
         ! L1 nodes i.e. nodes on L2 interpolated from L1, can only be hanging on L2 with the
         ! 1-irregularity index rule applied
         ! get nodes from all nodes with element level = meshLevel-1 i.e. these nodes are new
         ! nodes from L1 that form elements on L2
         ! write(*,*) "BaseMesh total vertices: ", base_mesh%numVertices
         do i = 1, base_mesh%numVertices
            vtx = base_mesh%v2pe(i)
            if ( vtx%elementLevel == parentLevel ) then
               ! this vertex is on `meshLevel`, check it
               ! get all cells incident on its parent edge in its elementLevel
               tets = base_mesh%getLevelElementsAttachedToEdge(parentLevel, vtx%parentEdge)
               ! if all are refined, not hanging node, if any is not refined, node is hanging
               ! write(*,*) "Size of edge attached level elements: ", size(tets)
               do j = 1, size(tets)
                  if ( .not. tets(j)%isNullElement() .and. tets(j)%isActive ) then
                     ! this tet is not refined
                     numHangingNodes = numHangingNodes + 1
                     tempNodes(i) = vtx
                     write(*,*) "hanging node is vtx with globalid: ", vtx%globalId
                     ! exit the loop cos hanging node found
                     exit
                  end if
               end do
            end if
            if ( vtx%elementLevel > parentLevel ) then
               exit
            end if
         end do
         ! check if hanging nodes were found
         if ( numHangingNodes > 0 ) then
            if ( allocated(hangingNodes) ) deallocate(hangingNodes)
            allocate(hangingNodes(numHangingNodes))
            ! extract these nodes from tempNodes
            ! to track nodes already slotted in hangingNodes
            k = 1
            do i = 1, base_mesh%numVertices
               vtx = tempNodes(i)
               if ( .not. vtx%isNullVertex() ) then
                  hangingNodes(k) = vtx
                  k = k + 1
               end if
            end do
            ! write(*,*) "Is K same as numHangingNodes: ", k==numHangingNodes
            ! write(*,*) "numHangingNodes: ", numHangingNodes
            ! write(*,*) "k: ", k
         end if
      end if
   end function getLevelMeshHangingNodes

   function createMidVertexOnEdge(base_mesh, edge, parentTet, vtxGlobalId) result(edgeMidVtx)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      type(MeshEdge), intent(in) :: edge
      type(Tetrahedron), intent(in) :: parentTet
      integer, intent(in) :: vtxGlobalId
      type(Vertex) :: edgeMidVtx, v1, v2
      double precision :: x, y, z

      v1 = base_mesh%v2pe(edge%vertices(1))
      v2 = base_mesh%v2pe(edge%vertices(2))
      x = (v1%x + v2%x) / 2
      y = (v1%y + v2%y) / 2
      z = (v1%z + v2%z) / 2

      edgeMidVtx = createVertex(x, y, z, vtxGlobalId, parentTet%elementLevelId, parentTet%level)
      ! assign `edge` as parentEdge for this mid vertex
      edgeMidVtx%parentEdge = edge
   end function createMidVertexOnEdge

   function getLevelMeshNumNodes(base_mesh, meshLevel) result(numLevelNodes)
      implicit none
      class(BaseMesh), intent(in) :: base_mesh
      integer, intent(in) :: meshLevel
      integer :: j, numLevelNodes
      type(Vertex) :: vtx

      ! to get numNodes for this level, count every node in order with parentLevel < meshLevel
      numLevelNodes = 0
      do j = 1, base_mesh%numVertices
         vtx = base_mesh%v2pe(j)
         if ( vtx%elementLevel < meshLevel ) then
            numLevelNodes = numLevelNodes + 1
         end if
         if(vtx%elementLevel >= meshLevel) exit
      end do
   end function getLevelMeshNumNodes

   !!!!! SUBROUTINES !!!!!!

   subroutine initialize(base_mesh)
      class(BaseMesh), intent(inout) :: base_mesh

      ! set stat values to zero
      base_mesh%numElements = 0
      base_mesh%numVertices = 0
      base_mesh%numLevels = 0

      ! allocate arrays
      ! let vertices start with space for 10
      allocate(base_mesh%vertices(10))
      allocate(base_mesh%v2pe(10))
      allocate(base_mesh%tetrahedrons(0))
      allocate(base_mesh%levelMeshes(0))
   end subroutine initialize

   subroutine initializeLevelMesh(base_mesh)
      class(BaseMesh), intent(inout) :: base_mesh
      type(LevelMesh) :: mesh
      integer :: meshLevel, i
      type(LevelMesh) :: levelMeshes(size(base_mesh%levelMeshes) + 1)

      base_mesh%numLevels = base_mesh%numLevels + 1
      levelMeshes(base_mesh%numLevels) = mesh
      do i = 1, size(base_mesh%levelMeshes)
         levelMeshes(i) = base_mesh%levelMeshes(i)
      end do
      base_mesh%levelMeshes = levelMeshes

      meshLevel = base_mesh%numLevels
      mesh = base_mesh%levelMeshes(meshLevel)
      mesh%numElements = 0
      mesh%level = meshLevel

      ! allocate(mesh%elementConn(0, 0))
      allocate(mesh%elements(0))
      allocate(mesh%sibhfcs(0,0))
      allocate(mesh%v2hf(0))
      allocate(mesh%v2es(0))
      allocate(mesh%e2pe(0))
      allocate(mesh%e2ce(0))

      base_mesh%levelMeshes(meshLevel) = mesh

   end subroutine initializeLevelMesh

   ! Add an array of tetrahedrons to the base Mesh and
   ! the levelMesh of the specified level
   ! A LevelMesh with level `meshLevel` MUST EXIST.
   ! Handles allocation of Cell Space
   subroutine addTetrahedrons(base_mesh, meshLevel, tetrahedrons)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      type(LevelMesh) :: mesh
      integer :: meshLevel, i
      ! type(Vertex) :: tetVertices(4)
      ! type(Vertex) :: hfVertices(3)
      type(Tetrahedron), allocatable :: tetrahedrons(:)


      ! allocate cells the size of the number of tetrahedrons in `tetrahedrons`
      call base_mesh%allocateCellSpace(meshLevel, size(tetrahedrons))
      mesh = base_mesh%levelMeshes(meshLevel)
      ! do for each tetrahedron in tetrahedrons
      do i = 1, size(tetrahedrons)
         base_mesh%numElements = base_mesh%numElements + 1
         mesh%numElements = mesh%numElements + 1
         ! add this vertices to base_mesh v2pe if mesh level > 1
         ! if ( meshLevel > 1 ) then
         !    ! New vertex from interpolation. map to parent
         !    call base_mesh%mapElementVerticesToParent(tetrahedrons(i))
         ! end if
         ! We're now adding original vertices from L1 to v2pe, giving them 0 as parent level Id and parent level
         ! this routine can run for every vertex now
         call base_mesh%mapElementVerticesToParent(tetrahedrons(i))
         ! add tet globalID to elements of level mesh
         ! tetrahedrons(i)%globalId = base_mesh%numElements
         mesh%elements(mesh%numElements) = tetrahedrons(i)%globalId

         ! add Tet to base_mesh tetrahedrons
         base_mesh%tetrahedrons(base_mesh%numElements) = tetrahedrons(i)
      end do

      ! save mesh to the base mesh
      base_mesh%levelMeshes(meshLevel) = mesh

   end subroutine addTetrahedrons

   subroutine expandV2pe(base_mesh)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      type(Vertex), allocatable :: v2pe(:)
      type(MeshConstants) :: meshConstant
      integer :: newSize, i
      ! expand it by 20% of its current size and copy the elements in it to the new one
      newSize = (size(base_mesh%v2pe) * meshConstant%v2peIncrement) + size(base_mesh%v2pe)

      allocate(v2pe(newSize))

      do i = 1, size(base_mesh%v2pe)
         v2pe(i) = base_mesh%v2pe(i)
      end do
      if ( allocated(base_mesh%v2pe) ) deallocate(base_mesh%v2pe)
      base_mesh%v2pe = v2pe
   end subroutine expandV2pe

   ! Handle new vertices, map vertex to ONE of its parents
   ! Add vertex to base Mesh. Expand v2pe. Basically saves the vertices
   ! As vertices are created, base_mesh%numVertices must keep incrementing to
   ! maintain the globalId count
   subroutine mapElementVerticesToParent(base_mesh, tet)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      type(Tetrahedron), intent(inout) :: tet
      integer :: i, globalId
      type(Vertex) :: vtx

      ! we want to assign global IDS to nodes not yet created here
      ! and attach to the first parent that sends it
      ! 4 nodes per Tet
      do i = 1, 4
         ! for each node, get the globalId
         globalId = tet%vertices(i)%globalId
         ! By checking if it's been saved already, we avoid updating the parentLevelId and the localId On current element
         ! if the globalId > size(v2pe) ==> it's never been saved before
         ! Automatically expand the v2pe or vertices array to twice its current size
         if ( globalId > size(base_mesh%v2pe) ) then
            call base_mesh%expandV2pe()
         end if

         vtx = base_mesh%v2pe(globalId)
         ! if v2pe(globalId)%globalId != globalId or == -1: this vertex has also not been saved before so:
         ! Update vertex localId on current element then save at v2pe(globalId) = Vertex<elementLevel, parentLevelId, localId>
         if ( vtx%globalId /= globalId .or. vtx%isNullVertex()) then
            tet%vertices(i)%localId = i
            base_mesh%v2pe(globalId) = tet%vertices(i)
            ! for every saved vertex, increment numVertices
            base_mesh%numVertices = base_mesh%numVertices + 1
         end if

      end do
   end subroutine mapElementVerticesToParent

   ! Adds tetrahedron `tet` to ALREADY RESERVED CELL SPACE in `base_mesh` and LevelMesh with level `meshLevel`
   subroutine addTetrahedron(base_mesh, meshLevel, tet)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      type(Tetrahedron), intent(inout) :: tet
      integer, intent(in) :: meshLevel
      type(LevelMesh) :: mesh
      ! integer :: j

      base_mesh%numElements = base_mesh%numElements + 1
      ! add this vertices to base_mesh v2pe if mesh level > 1
      if ( meshLevel > 1 ) then
         ! New vertex from interpolation. map to parent
         call base_mesh%mapElementVerticesToParent(tet)
      end if

      mesh = base_mesh%levelMeshes(meshLevel)
      mesh%numElements = mesh%numElements + 1

      ! add tet globalID to elements array of level mesh
      mesh%elements(mesh%numElements) = base_mesh%numElements
      ! add Tet to base_mesh tetrahedrons
      base_mesh%tetrahedrons(base_mesh%numElements) = tet

      ! save mesh to the base mesh
      base_mesh%levelMeshes(meshLevel) = mesh
   end subroutine addTetrahedron

   ! Expand the global array of tetrahedrons in base_mesh and the elementConn array on mesh with level `meshLevel` by `numCells`
   ! LEVELMESH MUST be initialized already, no error checking
   subroutine allocateCellSpace(base_mesh, meshLevel, numCells)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      integer :: meshLevel, numCells, i
      type(Tetrahedron), dimension(:), allocatable :: tetrahedrons
      ! integer, dimension(:,:), allocatable :: elementConn
      integer, allocatable :: elements(:)
      type(LevelMesh) :: mesh

      ! add numCells to their current sizes, allocate that and copy their previous entries
      allocate(tetrahedrons(size(base_mesh%tetrahedrons) + numCells))
      ! copy previous tetrahedrons
      do i = 1, size(base_mesh%tetrahedrons)
         tetrahedrons(i) = base_mesh%tetrahedrons(i)
      end do
      ! update base_mesh with new size
      if (allocated(base_mesh%tetrahedrons)) deallocate(base_mesh%tetrahedrons)
      base_mesh%tetrahedrons = tetrahedrons

      ! Get the Mesh at this level
      mesh = base_mesh%levelMeshes(meshLevel)
      ! 4 nodes form the connectivity
      allocate(elements(size(mesh%elements) + numCells))
      do i = 1, size(mesh%elements)
         elements(i) = mesh%elements(i)
      end do
      if (allocated(mesh%elements)) deallocate(mesh%elements)

      base_mesh%levelMeshes(meshLevel)%elements = elements

   end subroutine allocateCellSpace

   subroutine sortAdjacentVerticesOnFacet(vhfs)
      implicit none
      type(VertexHalfFacet), intent(inout) :: vhfs(2)
      type(VertexHalfFacet) :: vhf

      vhf = vhfs(1)
      ! We have just 2 adjacent vertices per face on a levelMesh
      ! if first vertex > 2nd vertex, swap both, otherwise do absolutely nothing
      if ( vhf%vertexGlobalId > vhfs(2)%vertexGlobalId  ) then
         vhfs(1) = vhfs(2)
         vhfs(2) = vhf
      end if
   end subroutine sortAdjacentVerticesOnFacet

   subroutine buildLevelMeshSibhfcs(base_mesh, meshLevel)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      type(LevelMesh) :: mesh
      integer :: meshLevel, i, j, k, vhfsCounter
      ! Intermediate Mappings to help get the sibling half-facets
      ! vertex with largest ID within facet => array of facets incident on the vertex
      ! in which vertex has largest ID. NO DUPLICATES OF VERTEX, half-facets are appended
      ! To retrieve vertex details in constant time, v2hfs will be the size of the
      ! number of vertices existing in base_mesh so,
      ! v2hfs(vertexGlobalID) = Vertex2HalfFacets struct with array of half-facets
      type(Vertex2HalfFacets), allocatable :: v2hfs(:)
      ! [
      !  v0 => [[VHF1_1, VHF1_2], [VHF2_1, VHF2_2]],
      !  v1 => [[V1HF1_1, V1HF1_2]]
      !]
      ! vertex (v) to its adjacent vertices in each of the half-facets in v2hfs(v)
      ! allocate size: base_mesh%numVertices for constant time retrieval with globalId
      type(Vertex2VertexHalfFacets), allocatable :: v2adj(:)
      ! To aid in internal v2hfs array expanding
      type(HalfFacet), allocatable :: hfs(:)
      ! To aid in internal v2adj array expanding
      ! array of vertex half-facets for set of adjacent vertices on a facet
      type(VertexHalfFacet) :: vhfs(2)
      type(Tetrahedron) :: tet
      type(HalfFacet) :: hf, sibhf
      type(Vertex) :: vtxMax
      type(VertexHalfFacet) :: vtxHf
      ! for each considered facet, hold vhfs(2) above => column size is always 2
      type(VertexHalfFacet), allocatable :: vhf(:,:)
      ! type(Vertex2HalfFacets) :: v2hf
      type(SiblingHFTuple) :: sibhfTuple
      type(SiblingHFTuple), allocatable :: tempSibhfcs(:,:)
      ! to aid with sorting
      ! integer, allocatable :: work(:)

      mesh = base_mesh%levelMeshes(meshLevel)

      ! sibhfcs actually needs to be updated if more elements are added per level
      ! we have to be careful to not entirely delete the data and reconstruct from scratch everytime
      ! we can do this by checking against the size
      write(*,*) "size of sibhfcs: ", size(mesh%sibhfcs, dim=1)
      if ( size(mesh%sibhfcs, dim=1) == 0) then
         ! allocate size of sibhfcs. num of elements in level mesh * num of facets per element
         if (allocated(mesh%sibhfcs)) deallocate(mesh%sibhfcs)
         allocate(mesh%sibhfcs(mesh%numElements, 4))
      end if

      ! check if the number of elements in the mesh has changed, if it has, sibhfcs needs to be
      ! expanded to accomodate the new size and the previous items copied
      if ( size(mesh%sibhfcs, dim=1) < mesh%numElements ) then
         allocate(tempSibhfcs(mesh%numElements, 4))
         ! copy current values in sibhfcs to the temp
         do i = 1, size(mesh%sibhfcs, dim=1)
            ! ! over the facets -> 4 facets
            ! do j = 1, 4
            !    tempSibhfcs
            ! end do
            tempSibhfcs(i,:) = mesh%sibhfcs(i,:)
         end do
         mesh%sibhfcs = tempSibhfcs
      end if

      allocate(v2hfs(base_mesh%numVertices))
      allocate(v2adj(base_mesh%numVertices))

      ! over the elements on this level
      do i = 1, mesh%numElements
         ! use tet's globalId to fetch already processed tet from
         tet = base_mesh%tetrahedrons(mesh%elements(i))
         ! over the half-facets in tet
         do j = 1, 4
            hf = tet%facets(j)
            ! get vertex with largest ID within this facet
            vtxMax = getHfVertexWithLargestId(hf)
            ! get set of adjacent vertices of `vtx` in this facet i.e. the remaining
            ! 2 vertices that are not this vtxMax
            ! over the vertices of this facet
            vhfsCounter = 1
            do k = 1, 3
               ! skip `vtxMax`
               if ( hf%vertices(k)%globalId == vtxMax%globalId ) cycle
               vtxHf%vertexGlobalId = hf%vertices(k)%globalId
               vtxHf%facetId = hf%localId
               vtxHf%elementGlobalId = tet%globalId
               vtxHf%elementLevelId = i
               vhfs(vhfsCounter) = vtxHf
               vhfsCounter = vhfsCounter + 1
            end do

            ! sort these adjacent vertices. This will make it easier to compare
            ! to other sibling half-facets as the ordering of the vertices will
            ! change in a sibling half-facet
            call sortAdjacentVerticesOnFacet(vhfs)

            ! Append facet into v2hfs(vtxMax)
            v2hfs(vtxMax%globalId)%vertexGlobalId = vtxMax%globalId
            ! allocate `hfs` array of `v2hfs(v)` if it hasn't been allocated before
            ! This ensures that FORTRAN doesn't give it a weird random integer
            ! to set it to zero. This will happen before the accumulation of hfs starts
            if (.not. allocated(v2hfs(vtxMax%globalId)%hfs)) allocate(v2hfs(vtxMax%globalId)%hfs(0))
            ! allocate hfs holder container. size of prev attached half-facets + 1
            ! hfs has to always be first deallocated before being allocated with new size
            if (allocated(hfs)) deallocate(hfs)
            allocate(hfs(size(v2hfs(vtxMax%globalId)%hfs) + 1))
            ! copy everything in the previous one if there's sth there i.e. the size is > 0
            if ( size(v2hfs(vtxMax%globalId)%hfs) > 0 ) then
               do k = 1, size(v2hfs(vtxMax%globalId)%hfs)
                  hfs(k) = v2hfs(vtxMax%globalId)%hfs(k)
               end do
            end if

            ! append this facet at last position
            hfs(size(v2hfs(vtxMax%globalId)%hfs) + 1) = hf
            if(allocated(v2hfs(vtxMax%globalId)%hfs))deallocate(v2hfs(vtxMax%globalId)%hfs)
            v2hfs(vtxMax%globalId)%hfs = hfs

            ! append set of adjacent vertices i.e. vhfs into v2adj(v,f)
            v2adj(vtxMax%globalId)%vertexGlobalId = vtxMax%globalId

            ! allocate the vhf array of v2adj to specify the numCols
            ! 2 cos there'd always be 2 adj vertices to the max vertice in every face
            if(.not. allocated(v2adj(vtxMax%globalId)%vhf)) allocate(v2adj(vtxMax%globalId)%vhf(0,2))
            ! allocate v2vhf holder container
            ! always 2 adjacent vertices per half-facet(face)
            ! write(*,*) "VHF row: ", size(v2adj(vtxMax%globalId)%vhf), "VHF col: ", size(v2adj(vtxMax%globalId)%vhf, dim=2)
            ! write(*,*) "VHF Data: ", v2adj(vtxMax%globalId)%vhf
            ! write(*,*) "v2adj size: ", size(v2adj)

            if (allocated(vhf)) deallocate(vhf)
            ! we need number of rows
            allocate(vhf(size(v2adj(vtxMax%globalId)%vhf, dim=1)+1, 2))
            ! copy previous entries
            do k = 1, size(v2adj(vtxMax%globalId)%vhf, dim=1)
               ! write(*,*) "vhf data for index: ", k
               ! write(*,*) v2adj(vtxMax%globalId)%vhf(k,:)
               vhf(k,:) = v2adj(vtxMax%globalId)%vhf(k,:)
            end do
            ! add entry last
            vhf(size(v2adj(vtxMax%globalId)%vhf, dim=1) + 1, :) = vhfs
            if(allocated(v2adj(vtxMax%globalId)%vhf))deallocate(v2adj(vtxMax%globalId)%vhf)
            v2adj(vtxMax%globalId)%vhf = vhf
         end do
      end do

      ! ***************************FOR DEBUGGING*****************************************
      ! write(*,*) "Currrent v2hfs: "
      ! do i = 1, size(v2hfs)
      !    write(*,*) "half-facets For Node ", i
      !    do k = 1, size(v2hfs(i)%hfs)
      !       write(*,*) "for half-facet ", k
      !       write(*,*) "facet element level ID: ", v2hfs(i)%hfs(k)%elementLevelId
      !       write(*,*) "facet local Id: ", v2hfs(i)%hfs(k)%localId
      !       write(*,*) "facet's vertices: ", v2hfs(i)%hfs(k)%vertices
      !    end do
      ! end do

      ! write(*,*) "Currrent V2adj: "
      ! do i = 1, size(v2adj)
      !    write(*,*) "Adjacent vertices in half-facets for Node ", i
      !    ! over the half-facets incident on the vertice
      !    do k = 1, size(v2adj(i)%vhf, dim=1)
      !       write(*,*) "for half-facet ", k
      !       ! over the exactly 2 adjacent vertex-halffacet pair per half-facet
      !       do j = 1, 2
      !          write(*,*) "Adjacent vertex-Hf pair ", j
      !          write(*,*) "Confirm vertex global ID: ", v2adj(i)%vertexGlobalId
      !          write(*,*) "Adjacent Vertex global ID", v2adj(i)%vhf(k,j)%vertexGlobalId
      !          write(*,*) "Attached facet local Id: ", v2adj(i)%vhf(k,j)%facetId
      !          write(*,*) "Facet element's level Id: ", v2adj(i)%vhf(k,j)%elementLevelId
      !          write(*,*) "Facet element's global Id: ", v2adj(i)%vhf(k,j)%elementGlobalId
      !       end do

      !    end do
      ! end do
      ! Again over all elements in level mesh. for actual sibhfcs
      ! sibhfcs for level mesh should be already set here, so there'd be no need
      ! to rebuild it entirely after the first time it's been built for this levelMesh
      do i = 1, mesh%numElements
         ! use tet's globalId to fetch already processed tet from
         tet = base_mesh%tetrahedrons(mesh%elements(i))
         ! over the facets of this element
         do j = 1, 4
            hf = tet%facets(j)
            ! sibhfcs(element, facet)
            if ( .not. isSiblingHFTupleSet(mesh%sibhfcs(i,j)) ) then
               vtxMax = getHfVertexWithLargestId(hf)
               ! get set of adj vertices of `vtxMax` in v2adj(v,f)
               vhfs = getAdjacentVerticesToMaxVtxInFacet(v2adj, vtxMax%globalId, hf)
               ! find half-facets in v2hfs(vtxMax) such that v2adj(v,hfi)=vhfs
               ! There should just be a max of 2 matching facets, the current considered
               ! facet and its sibling, would be 1 if this facet has no sibling then just it will be found

               hfs = getHalfFacetsWithSameAdjacentVertices(v2hfs, vtxMax%globalId, v2adj, vhfs)
               ! Force hfs to return the number of facets found, either one or two
               ! use the size to check if there's a sibling. 2 means sibling found, 1, no sibling
               ! Form cyclic mapping btw the returned half-facets in sibhfs
               ! i.e. hf1 <--> hf2
               if ( size(hfs) > 1 ) then
                  ! should have a sibling
                  sibhf = hfs(1)
                  if ( hf%isHfEqual(sibhf) ) then
                     sibhf = hfs(2)
                  end if
                  sibhfTuple%sibhfElemLevelId = sibhf%elementLevelId
                  sibhfTuple%sibhfLocalId = sibhf%localId
                  ! setting this facet to its sibling
                  mesh%sibhfcs(i,j) = sibhfTuple
                  ! setting for the sibling facet to this facet
                  sibhfTuple%sibhfElemLevelId = hf%elementLevelId
                  sibhfTuple%sibhfLocalId = hf%localId
                  mesh%sibhfcs(sibhf%elementLevelId, sibhf%localId) = sibhfTuple
               else
                  ! should have no sibling. Set its Sibling to NULL sibhfTuple
                  sibhfTuple%sibhfElemLevelId = 0
                  sibhfTuple%sibhfLocalId = 0
                  mesh%sibhfcs(i,j) = sibhfTuple
               end if

            end if
         end do
      end do

      ! update the base_mesh with new updates to the level mesh
      base_mesh%levelMeshes(meshLevel) = mesh

   end subroutine buildLevelMeshSibhfcs

   ! To construct mapping from each vertex to an incident half-facet
   subroutine buildLevelMeshV2hf(base_mesh, meshLevel)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      integer, intent(in) :: meshLevel
      type(HalfFacet) :: hf
      integer :: i, j, k
      type(LevelMesh) :: mesh
      type(Tetrahedron) :: tet
      type(Vertex) :: vtx

      mesh = base_mesh%levelMeshes(meshLevel)
      ! allocate level mesh v2hf, for 0(1) access, we're using an array as a HashMap
      ! We allocate v2hf to the size of nodes in base_mesh, so the index of v2hf will
      ! be the globalId of the vertex. v2hf(v%globalId) => HalfFacet
      if(allocated(mesh%v2hf)) deallocate(mesh%v2hf)
      allocate(mesh%v2hf(base_mesh%numVertices))

      ! over the elements
      do i = 1, mesh%numElements
         ! We store the global IDs of the elements in the levelMesh elements array
         ! to retrieve from the global base_mesh
         tet = base_mesh%tetrahedrons(mesh%elements(i))
         ! over the vertices in Element
         do j = 1, 4
            vtx = tet%vertices(j)
            ! map v2hf(v) to a facet incident on v in tet if v2hf(v) is not set
            ! To do this: recall the CGNS labelling, from this, v1 can never be on facet3,
            ! v2 can never be on facet4, v3 never on facet2, v4 never on facet1
            ! Pattern is v_i is always on facet_i
            if ( mesh%v2hf(vtx%globalId)%isNullHf() ) then
               mesh%v2hf(vtx%globalId) = tet%facets(j)
            end if
         end do
         ! To prioritize border facets. Border facets are such that sibhfcs<tet,hf> is NULL, no sibling
         ! Border facets are prioritized to be able to use v2hf(v) to determine if v is a border vertex
         ! if facet is on the border => vertex on it is also a border vertex
         ! over the facets in tet
         do j = 1, 4
            hf = tet%facets(j)
            ! if facet is a border facet, map its vertices to itself
            if ( mesh%sibhfcs(i, j)%isNullSibHFTuple() ) then
               ! over the vertices of hf
               do k = 1, 3
                  vtx = hf%vertices(k)
                  mesh%v2hf(vtx%globalId) = hf
               end do
            end if
         end do
      end do
      ! update the base_mesh with new updates to the level mesh
      base_mesh%levelMeshes(meshLevel) = mesh
   end subroutine buildLevelMeshV2hf

   subroutine buildLevelMeshV2es(base_mesh, meshLevel)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      integer, intent(in) :: meshLevel
      type(LevelMesh) :: mesh
      ! HashSet to accumulate element global IDs incidented on some vertex
      ! type(ftlHashSetInt) :: tempSet
      integer, allocatable :: tempArr(:)
      integer :: i, j, k
      type(Tetrahedron) :: tet
      type(Vertex) :: vtx
      type(Vertex2Elements) :: v2elem

      mesh = base_mesh%levelMeshes(meshLevel)

      ! we clear v2es and rebuild it from scratch each time cos we want to avoid
      ! duplicate element values and the waste of time to do a linear search of an element
      ! in the array. This can simply be solved by using a hashset but fortran is limited
      write(*,*) "size of v2es: ", size(mesh%v2es)
      ! allocate size of v2es. num of nodes in base_mesh
      ! using an array as a map.
      ! v2es(v%globalId) = set of integers(localIds of elements on this level)
      if (allocated(mesh%v2es)) deallocate(mesh%v2es)
      allocate(mesh%v2es(base_mesh%numVertices))
      ! for every element on this level, for every vertex in this element
      ! allocate size of tempArr to size(mesh%v2es(vtx%globalId)+1), copy into tempArr
      ! deallocate mesh%v2es(vtx%globalId), set to tempArr
      do i = 1, mesh%numElements
         ! Global Ids of tetrahedrons stored per level to access it from base_mesh
         tet = base_mesh%tetrahedrons(mesh%elements(i))
         do j = 1, 4
            vtx = tet%vertices(j)
            v2elem = mesh%v2es(vtx%globalId)
            v2elem%vertexGlobalId = vtx%globalId
            ! if v2elems array have never had a prior allocation, allocate size 0
            if (.not. allocated(v2elem%v2elems)) allocate(v2elem%v2elems(0))
            if (allocated(tempArr)) deallocate(tempArr)
            allocate(tempArr(size(v2elem%v2elems) + 1))
            ! copy to tempArr
            do k = 1, size(v2elem%v2elems)
               tempArr(k) = v2elem%v2elems(k)
            end do
            ! add new element last
            tempArr(size(v2elem%v2elems) + 1) = tet%elementLevelId
            if (allocated(v2elem%v2elems)) deallocate(v2elem%v2elems)
            v2elem%v2elems = tempArr
            ! update
            mesh%v2es(vtx%globalId) = v2elem
         end do
      end do
      ! update base_mesh
      base_mesh%levelMeshes(meshLevel) = mesh
   end subroutine buildLevelMeshV2es

   subroutine buildLevelMeshInterLevelConnectivity(base_mesh, meshLevel)
      class(BaseMesh), intent(inout) :: base_mesh
      integer, intent(in) :: meshLevel
      integer, allocatable :: tempArr(:)
      type(LevelMesh) :: mesh
      integer :: i
      type(Tetrahedron) :: tet

      mesh = base_mesh%levelMeshes(meshLevel)
      ! To update, e2pe and e2ce or just keep their size same as numElements on level mesh

      ! for every element on the level, get the parentElementId and slot in e2pe if level >1
      ! level 1 don't have parents, we don't want to waste space storing 0s or -1s
      allocate(tempArr(mesh%numElements))
      if ( meshLevel > 1 ) then
         do i = 1, mesh%numElements
            tet = base_mesh%tetrahedrons(mesh%elements(i))
            tempArr(i) = tet%parentElementLevelId
         end do
         if ( allocated(mesh%e2pe) ) deallocate(mesh%e2pe)
         mesh%e2pe = tempArr
      end if
      if ( size(mesh%e2ce) /= mesh%numElements ) then
         ! copy into tempArr
         if(allocated(tempArr)) deallocate(tempArr)
         allocate(tempArr(mesh%numElements))

         ! No consistent default value for ints in fortran and e2ce has to be checked
         ! No need to stress by directly checking the e2ce, if element is inactive i.e.
         ! it's been refined, the e2ce(elem) will be correct

         do i = 1, size(mesh%e2ce)
            tempArr(i) = mesh%e2ce(i)
         end do
         if(allocated(mesh%e2ce)) deallocate(mesh%e2ce)
         mesh%e2ce = tempArr
      end if
      base_mesh%levelMeshes(meshLevel) = mesh
   end subroutine buildLevelMeshInterLevelConnectivity

   subroutine buildLevelMeshConnectivity(base_mesh, meshLevel)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      integer, intent(in) :: meshLevel
      ! type(LevelMesh) :: mesh

      ! mesh = base_mesh%levelMeshes(meshLevel)

      ! build sibling half-facets mapping for the level Mesh for adjacency info
      call base_mesh%buildLevelMeshSibhfcs(meshLevel)
      ! base_mesh%levelMeshes(meshLevel) = mesh
      ! build optional vertex to half-facet mapping for level Mesh
      ! This is needed to easily check if a node is a border node on a level
      call base_mesh%buildLevelMeshV2hf(meshLevel)

      ! vertex to elements incident on it on level `meshLevel`
      call base_mesh%buildLevelMeshV2es(meshLevel)

      ! e2pe & e2ce need to be kept updated, to be the size of the numElements on that LevelMEsh
      call base_mesh%buildLevelMeshInterLevelConnectivity(meshLevel)
   end subroutine buildLevelMeshConnectivity

   ! refine element with the 1-irregularity rule
   ! will recursively refine parent neighbours if they're unrefined to achieve a graded mesh
   recursive subroutine refineLevelElement(base_mesh, meshLevel, elementLevelId)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      integer, intent(in) :: meshLevel, elementLevelId
      type(LevelMesh) :: mesh, parentMesh
      integer :: i, j, parentLevel, midVtxCount, numNodes, nextLevel
      type(Tetrahedron) :: tet, refinedTet
      type(Tetrahedron), allocatable :: largeNeighbours(:), edgeNeighbours(:)
      type(Vertex) :: midVertices(6)
      type(Vertex) :: edgeMidVtx

      write(*,*) "Attempting to refine element with level ID: ", elementLevelId, "on level: ", meshLevel
      mesh = base_mesh%levelMeshes(meshLevel)
      tet = base_mesh%tetrahedrons(mesh%elements(elementLevelId))
      parentLevel = meshLevel - 1
      nextLevel = meshLevel + 1
      ! exit if element is already refined. Base case for the recursive subroutine
      if ( .not. tet%isActive ) then
         write(*,*) "Element already refined!"
         return
      end if

      ! Confirm that the children level exists in base_mesh, if it doesn't
      ! create and initialize it
      if ( base_mesh%numLevels < nextLevel ) then
         call base_mesh%initializeLevelMesh()
      end if

      ! if element's meshLevel > 1 then check its large neigbours to ensure they've been
      ! refined. These are elements on its parent level that it shares edges with.
      ! level 1 elements are original, no parents
      if ( meshLevel > 1 ) then
         parentMesh = base_mesh%levelMeshes(parentLevel)
         ! get Large neighbours for tet
         largeNeighbours = base_mesh%getElementLargeNeighbours(tet)
         ! over the returned neighbours
         do i = 1, size(largeNeighbours)
            ! No NULL element returned here, just check if it's unrefined
            if (largeNeighbours(i)%isActive) then
               write(*,*) "Unrefined large neighbour found, with details: "
               write(*,*) "Element Level Id: ", largeNeighbours(i)%elementLevelId
               write(*,*) "On Level: ", parentLevel
               call base_mesh%refineLevelElement(parentLevel, largeNeighbours(i)%elementLevelId)
            end if
         end do
      end if

      ! At this point, every facet-neighbour of this element should be refined
      ! now, to either use existing midpoint nodes or create new ones if they dont exist
      ! Loop over the edges in 'tet' to refine the edge or extract midpoint from already
      ! refined edge

      midVtxCount = 1
      numNodes = base_mesh%numVertices
      ! over tet's edges
      do i = 1, 6
         write(*,*) "Refining edge ", i
         ! get every cell incidented on the edge
         edgeNeighbours = base_mesh%getLevelElementsAttachedToEdge(tet%level, tet%edges(i))
         ! Loop over the cells to find one that is refined to then access it for its midedge node
         do j = 1, size(edgeNeighbours)
            refinedTet = edgeNeighbours(j)
            ! not active => it's already refined
            if ( .not. refinedTet%isActive ) exit
         end do
         ! see CGNS system, the midpoints are extracted/created in ORDER. vtx 5 - 10
         ! check this already refined neighbour that shares the same edge for the mid node
         edgeMidVtx = base_mesh%getEdgeMidVertex(tet%edges(i), refinedTet)
         if ( edgeMidVtx%isNullVertex() ) then
            write(*,*) "No already existing mid vertex, creating mid vertex..."
            ! mid vtx not found, create it
            numNodes = numNodes + 1
            edgeMidVtx = base_mesh%createMidVertexOnEdge(tet%edges(i), tet, numNodes)
            write(*,*) "Parent edge of mid vtx: ", edgeMidVtx%parentEdge%vertices
         end if
         write(*,*) "Edge v1->v2: ", tet%edges(i)%vertices
         write(*,*) "Mid-vertex for above edge has globalId: ", edgeMidVtx%globalId
         midVertices(midVtxCount) = edgeMidVtx
         midVtxCount = midVtxCount + 1
      end do
      ! base_mesh%numVertices will eventually be updated after the refined children of tet
      ! are added to base_mesh
      ! at this point, all midpoint nodes should have been gotten and duplicates avoided
      ! so, split
      call base_mesh%localRegularRefinement(tet, midVertices)
      ! build mesh connectivity for children level
      call base_mesh%buildLevelMeshConnectivity(meshLevel + 1)
   end subroutine refineLevelElement

   subroutine localRegularRefinement(base_mesh, tet, midVertices)
      implicit none
      class(BaseMesh), intent(inout) :: base_mesh
      type(Tetrahedron), intent(in) :: tet
      type(Vertex), intent(in) :: midVertices(6)
      type(Tetrahedron), allocatable :: subTets(:)
      integer :: numLevelElements, numElements, nextLevel, parentLevelId
      type(Vertex) :: vertices(4)

      write(*,*) "Splitting tet with element level ID: ", tet%elementLevelId, "on level: ", tet%level

      allocate(subTets(8))
      nextLevel = tet%level + 1
      parentLevelId = tet%elementLevelId

      numElements = base_mesh%numElements
      numLevelElements = base_mesh%levelMeshes(nextLevel)%numElements

      ! J.Bey's refinement algorithm
      ! midEdgeNodes 5 to 10 are midVertices 1 to 6
      ! T1 = 1,5,7,8
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/tet%vertices(1), midVertices(1), midVertices(3), midVertices(4)/)
      subTets(1) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T2 = 5,2,6,9
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(1), tet%vertices(2), midVertices(2), midVertices(5)/)
      subTets(2) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T3 = 7,6,3,10
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(3), midVertices(2), tet%vertices(3), midVertices(6)/)
      subTets(3) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T4 = 8,9,10,4
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(4), midVertices(5), midVertices(6), tet%vertices(4)/)
      subTets(4) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T5 = 5,7,8,9
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(1), midVertices(3), midVertices(4), midVertices(5)/)
      subTets(5) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T6 = 5, 7, 6, 9
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(1), midVertices(3), midVertices(2), midVertices(5)/)
      subTets(6) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T7 = 7,8,9,10
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(3), midVertices(4), midVertices(5), midVertices(6)/)
      subTets(7) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      ! T8 = 7, 6, 9, 10
      numElements = numElements+1
      numLevelElements = numLevelElements + 1
      vertices = (/midVertices(3), midVertices(2), midVertices(5), midVertices(6)/)
      subTets(8) = createTetrahedron(vertices, numLevelElements, numElements, nextLevel, parentLevelId)

      write(*,*) "Adding the child tets to next level of base_mesh..."
      call base_mesh%addTetrahedrons(nextLevel, subTets)
      write(*,*) "Child Tets added successfully!"
      ! update e2ce for tet here
      write(*,*) "Updating e2ce for parent element..."
      ! first child element level Id
      base_mesh%levelMeshes(tet%level)%e2ce(tet%elementLevelId) = subTets(1)%elementLevelId
      ! include parentLevelId for subTets

      ! update parent isActive flag to false cos it's been refined
      base_mesh%tetrahedrons(tet%globalId)%isActive = .false.
   end subroutine localRegularRefinement

end module MeshLib
