# Information about the MUI Linear Algebra
This is MUI linear algebra system with namespace of mui::linalg.
The MUI linear algebra system support COO, CSR and CSC matrix formats.
Linear solvers and preconditioners have been implemented.

# Compilation notes
It relies on C++11 standard.

# Folder and file structure
The MUI linear algebra system with namespace of mui::linalg is self-contained. Three top level files are in the folder on different aspects:
* matris.h
-- It is the base file for the MUI sparse matrix system.
* solver.h
-- It is the base file of the linear system solver.
* preconditioner.h
-- It is the based file of the preconditioners for the linear system solver.
* linalg_util.h
-- It is the file of the utility functions for the MUI linear algebra system.

Implementations are in seperate files with the name of the basic file name plus underscore and the name of the specific implementation content.

# Potential performance improvement
- Current MUI linear system solver implemented Gaussian Elimination (direct) solver, Conjugate Gradient (iteritive) solver and BiCGStab solver. Generalized minimal residual method (GRMRES) with restart function should be implemented for large matrix.

- Current MUI linear algebra system are mainly serial codes, which makes the assumption that MUI has already finished the particion of the matrix (points) into a number of processors. The MUI linear algebra system, then, be applied into every processor to mainipulate/solver the linear formulations locally. If there are scenarios that needed other processors to help the linear formulations mainipulation/solving of the local processor, the serial MUI linear algebra code should be extended to parallel. 

# Unit tests
Unit test code are provided in the "test" folder:
1. matrix_arithmetics.cpp and matrix_manipulations.cpp
They are test codes for the MUI sparse matrix system, covers the test of basic matrix setting, I/O, mainipulation, arithmetic operations (addition, subtraction, multiplication, scalar multiplication, hadamard product, transpose and inverse), LU decomposition and QR decomposition.
2. solver.cpp
It is a test code for the MUI linear system solvers, covers the MUI Gaussian Elimination solver and MUI Conjugate Gradient solver with different preconditioners (Incomplete LU, Incomplete Cholesky and Symmetric Successive Over Relaxation) and with/without initial guess.