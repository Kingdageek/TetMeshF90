# TetMeshF90

Fortran mesh Library, based on the Hierarchical Array-based Half-facet Data Structure by X. Zhao et al, for representation and multilevel refinement of 3D tetrahedral meshes.

This library was developed as part of my thesis for an M.Sc. in Material Science & Engineering at the University of Kiel titled _Testing and Extension of Multigrid Finite Element Method (FEM) Code: Implementation of volume meshing with tetrahedrons and mesh refinement_. As such, it's still experimental software. You can read the thesis here: [https://drive.google.com/file/d/1QCmXw6_mXAnN-ntSbp22KcJz1R90CftN/view?usp=sharing]. You have to read at least chapter 3 on Methodology.
The underlying mesh data structure is based on the work published in these two papers:

- V. Dyedov, N. Ray, D. Einstein, and X. Jiao, “AHF : Array-based half-facet data structures for mixed-dimensional and non- manifold meshes . AHF : Array-based Half-Facet Data Structure for Mixed-Dimensional and Non-manifold Meshes,” no. May 2014, 2013.
- X. Zhao, R. Conley, N. Ray, V. S. Mahadevan, and X. Jiao, “Conformal and non-conformal adaptive mesh refinement with hierarchical array-based half-facet data structures,” Procedia Eng., vol. 124, pp. 304–316, 2015, doi: 10.1016/j.proeng.2015.10.141.

You don't have to read them though, unless you want to, as chapter 3 in my thesis already describes in simple terms the conventions used and the underlying data structure for the library.

The library is modular and provides functions and subroutines that can be used to perform Adaptive Mesh Refinement for tetrahedral meshes.

### Current Features

The library can be currently used to:

- Represent an arbitrary conformal mesh of tetrahedrons
- Perform Uniform and Adaptive Mesh Refinement (Each Tetrahedron is regularly refined - 8 sub-tetrahedrons produced - red refinement based on J. Bey's algorithm)
- Track elements, nodes, and hanging nodes per level
- Generate files for use in computation say for a FEM procedure: can generate elements file per level, nodes file, hanging nodes file.
- Generate `.vtu` files for mesh visualization using ParaView
- Locate boundary nodes
- Retrieve adjacency information: elements directly abutting the faces of an element can be retrieved, every element around an element can be retrieved, every element around an edge can be retrieved, and lastly every large neighbour element one level higher than an element's level can also be retrieved.

### Restrictions

These restrictions are imposed:

- The initial, arbitrary, coarse mesh must be a conformal and manifold 3D mesh of tetrahedrons
- There cannot be more than one hanging node per shared edge (or face) of an element. This is the 1-irregularity index rule. This is automatically enforced during mesh refinement and allows for the generation of graded meshes

### Use

You have to create your own FORTRAN programs for what you want to implement and then use the procedures in the library to represent the meshes and refine them. Also, compile the library into your program to run. The appendix section in the thesis document details how to do all this. Also, see the _src/regular_tet.f90_ file [https://github.com/Kingdageek/TetMeshF90/blob/master/src/regular_tet.f90] for a sample code on how it's used.

Also, take a look at the _BaseMesh_ structure in the _src/MeshLib.f90_ file [https://github.com/Kingdageek/TetMeshF90/blob/master/src/MeshLib.f90], for an overview of all the procedures available.
