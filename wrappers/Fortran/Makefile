# Define compilers and MPI run option
FC=mpifort
CCX=mpic++
MPI=mpirun

# C compilation flags
FFLAGS = -O3 -cpp -std=f2003 -fcheck=all -Wall -g
# C++ compilation flags
CXXFLAGS = -O3 -std=c++11 -fpic -g

ifdef USE_RBF
# Path to Eigen library
EIGEN_PATH = /home/stephen/eclipse-workspace/eigen-3.3.9
# C compilation flags
FFLAGS += -DUSE_RBF
# C++ compilation flags
CXXFLAGS += -DUSE_RBF -I${EIGEN_PATH}
endif

default: mui_f_wrapper_1d.o
	 @echo "Compiling unit test code..."
	 ${FC} ${FFLAGS} mui_f_wrapper_1d.f90 unit_test.f90 mui_f_wrapper_1d.o -o unit_test_fortran_wrapper.x -lstdc++ -lmpi_cxx
	 @echo "Launching unit test code..."
	 ${MPI} -np 1 ./unit_test_fortran_wrapper.x domain1 interface : -np 1 ./unit_test_fortran_wrapper.x domain2 interface
	 @echo "************************************************************"
	 @echo "*              Fortran wrapper usage                       *"
	 @echo "*                                                          *"
	 @echo "* Step 1) Include compiled 1D/2D/3D module                 *"
	 @echo "* Step 2) Compile with f2003 compliant compiler            *"
	 @echo "************************************************************"

mui_f_wrapper_1d.o: mui_f_wrapper_1d.cpp
	 @echo "Generating Fortran-wrapper 1D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_f_wrapper_1d.cpp -o mui_f_wrapper_1d.o

clean:
	 rm -rf *.o *.mod unit_test_fortran_wrapper.x
