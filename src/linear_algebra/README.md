# Information about the MUI Linear Algebra
This is MUI linear algebra system with namespace of mui::linalg.

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
- The current sparse matrix system uses std::map(std::pair(ITYPE,ITYPE),VTYPE) to store matrix non-zero elements, where std::pair(ITYPE,ITYPE) stores the element index (i.e. row and column number) and VTYPE stores the value of the element. std::map() is implemented as a red-black tree, which provides logarithmic time complexity for both inserting and accessing elements, making it suitable for sparse matrices with a small number of non-zero elements. However, std::map also has additional overhead for managing the tree structure, and the cost of searching for an element can be high. Potential replace std::map() into std::vector, which provides constant time complexity for accessing elements. When the number of non-zero elements is large and the indices are in a contiguous range, std::vector can be more efficient.

- The current sparse matrix system uses COO format to stroe the matrix. Potential replace COO format into CSR or CSC formats for better performance.

- Current MUI linear system solver implemented Gaussian Elimination (direct) solver and Conjugate Gradient (iteritive) solver. Generalized minimal residual method (GRMRES) with restart function should be implemented for large matrix.

- Current MUI linear algebra system are mainly serial codes, which makes the assumption that MUI has already finished the particion of the matrix (points) into a number of processors. The MUI linear algebra system, then, be applied into every processor to mainipulate/solver the linear formulations locally. If there are scenarios that needed other processors to help the linear formulations mainipulation/solving of the local processor, the serial MUI linear algebra code should be extended to parallel. 

# Unit tests
Unit test code are provided in the "test" folder:
1. matrix.cpp
It is a test code for the MUI sparse matrix system, covers the test of basic matrix setting, I/O, mainipulation, arithmetic operations (addition, subtraction, multiplication, scalar multiplication, hadamard product, transpose and inverse), LU decomposition and QR decomposition.
2. solver.cpp
It is a test code for the MUI linear system solvers, covers the MUI Gaussian Elimination solver and MUI Conjugate Gradient solver with different preconditioners (Incomplete LU, Incomplete Cholesky and Symmetric Successive Over Relaxation) and with/without initial guess.