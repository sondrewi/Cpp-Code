# Cpp-Code
A selection of Cpp files and associated text files created as part of the course Numerical Methods for Incompresible Fluid Dynamics (Scientific Computing MPhil, Cambridge University).
Given the Mesh class, the sparse matrix addressing class (SparseAddress.cpp) constructs an object which stores the matrix structure of a mesh object for discretisation of a differential equation over the associated domain.
This is done assuming a compressed row format. 
The SparseMat.cpp class creates a sparse matrix object derived from a SparseAddress object. When populated with suitable coefficients, the SparseMat class can determine various attributes of the sparse matrix, such as symmetry,
positive definiteness, and determinants of submatrices. 
The Sparse Linear System class (SLS.cpp) creates a sparse linear system derived from a sparse matrix object, A. Furthermore, it stores vectors x and b appearing in the linear system Ax = b. The SLS.cpp-file contains various
functions used to actually discretise the terms of the convection-diffusion equation and to thereby populate the entries of the sparse matrix object A and the vector b.
The Gauss-Seidel Smoother class (GSSmooth.cpp) defines an object derived from a sparse linear system object. A GSSmooth object has a function "smooth" which runs through the provided number of Gauss-Seidel 
sweeps in order to solve the system Ax = b as defined in the sparse linear system object from which it was derived. 
Finally, the Gauss-Seidel Solver class (GSSolver.cpp) defines objects derived from a GSSmooth object. A GSSolver object has a function "Solve", which will call on the "smooth" function from its GSSmooth object until
the residual of the associated sparse linear system is sufficiently small or until a maximum number of iterations has been reached.

The GSSolverTest.cpp file defines a mesh (whose source files are contained in mesh_files) and associated sparse address, sparse matrix, sparse system, smoother, and solver objects in order to solve the convection-diffusion
equation over the mesh domain. For the supplied source files this is a rectangular 2D 3m*0.5m domain. At an initial time zero, GSSolverTest.cpp uniformly specifies some scalar T to be 1 throughout the entire domain.
A left boundary Dirichlet condition of T=0 and Neumann boundary conditions grad(T)=0 are also defined alongside an intermittent source term at the location (0.5,0.0). GSSolverTest.cpp then solves the convection-diffusion
equation for 1000 timesteps (with each timestep being 0.01 second). The computed temporal evolution of scalar T is then shown in the animation contained in ConvectionDiffusion.gif. 
