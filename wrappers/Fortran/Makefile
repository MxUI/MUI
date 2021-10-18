# Define compilers and MPI run option
FC=mpifort
CCX=mpic++
MPI=mpirun
AR=ar

# C compilation flags
FFLAGS = -O3 -std=f2003 -cpp -g -fcheck=all -Wall
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
	 ${FC} ${FFLAGS} mui_f_wrapper_1d.F90 unit_test.F90 mui_f_wrapper_1d.o -o unit_test.x -lstdc++ -lmpi_cxx
	 @echo "************************************************************"
	 @echo "*              Fortran wrapper usage                       *"
	 @echo "*                                                          *"
	 @echo "* Step 1) include mui_c_wapper_[1d/2d/3d].h in C source    *"
	 @echo "* Step 2) Compile with a C compiler & link against library *"
	 @echo "************************************************************"
	 
mui_f_wrapper_1d.o: mui_f_wrapper_1d.cpp
	 @echo "Generating Fortran-wrapper 1D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_f_wrapper_1d.cpp -o mui_f_wrapper_1d.o

clean:
	 rm -rf *.o *.mod unit_test_fortran_wrapper.x
