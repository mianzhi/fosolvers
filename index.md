**`foeuler`, A Transient Euler Equations (Inviscid Flow) Solver**

Here is a very brief [user guide](https://drive.google.com/file/d/0B0ArQPxS-HQZWktBSzIwSlRzdk0/view?usp=sharing)

* 3-D unstructured grid
* [GMSH](http://geuz.org/gmsh/) for grid generation
* [ParaView](http://www.paraview.org/) for post-processing
* Solves [internal](https://drive.google.com/file/d/0B0ArQPxS-HQZOFE3NmxWd296Vlk/view?usp=sharing) and [external](https://drive.google.com/file/d/0B0ArQPxS-HQZdWJQcEdxRENMX2M/view?usp=sharing) flow problems
* 1-st order Roe's scheme
* BDF-Newton-GMRES [solver](https://computation.llnl.gov/casc/sundials/main.html)
* Threaded Jacobi iteration preconditioner
* Simple UDF by [libmatheval](http://www.gnu.org/software/libmatheval/)

***

**`libfosolvers`, A Library for Scientific Computation Programs**

* Data structures and procedures for polyhedron grid
* Data structures and procedures for octree grid
* Geometric computation for polyhedron
* Grid input, data output procedures
* Finite volume spatial schemes
* General purpose condition management framework
* General purpose UDF framework by [libmatheval](http://www.gnu.org/software/libmatheval/)