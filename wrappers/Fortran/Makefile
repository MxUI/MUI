# Define compilers and MPI run option
FC=mpifort
CCX=mpic++
MPI=mpirun
AR=ar

# C compilation flags
FFLAGS = -O3 -cpp -std=f2003 -fcheck=all -Wall
# C++ compilation flags
CXXFLAGS = -O3 -std=c++11 -fpic -shared

ifdef USE_RBF
    # Path to Eigen library
    EIGEN_PATH = /usr/include/eigen3
    # C compilation flags
    FFLAGS += -DUSE_RBF
    # C++ compilation flags
    CXXFLAGS += -DUSE_RBF -I${EIGEN_PATH}
endif

default: mui_f_wrapper_general.o mui_f_wrapper_1d.o mui_f_wrapper_2d.o mui_f_wrapper_3d.o mui_f_wrapper_general.mod mui_f_wrapper_1d.mod mui_f_wrapper_2d.mod mui_f_wrapper_3d.mod
	 @echo "Compiling and linking MUI Fortran wrapper..."
	 ${CCX} ${CXXFLAGS} -o libmui_f_wrapper.so mui_f_wrapper_general.o mui_f_wrapper_1d.o mui_f_wrapper_2d.o mui_f_wrapper_3d.o
	 ${AR} rcs libmui_f_wrapper.a mui_f_wrapper_general.o mui_f_wrapper_1d.o mui_f_wrapper_2d.o mui_f_wrapper_3d.o

mui_f_wrapper_general.o: mui_f_wrapper_general.cpp
	 @echo "Generating MUI Fortran wrapper general object file..."
	 ${CCX} ${CXXFLAGS} -c mui_f_wrapper_general.cpp -o mui_f_wrapper_general.o

mui_f_wrapper_1d.o: mui_f_wrapper_1d.cpp
	 @echo "Generating MUI Fortran wrapper 1D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_f_wrapper_1d.cpp -o mui_f_wrapper_1d.o

mui_f_wrapper_2d.o: mui_f_wrapper_2d.cpp
	 @echo "Generating MUI Fortran wrapper 2D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_f_wrapper_2d.cpp -o mui_f_wrapper_2d.o

mui_f_wrapper_3d.o: mui_f_wrapper_3d.cpp
	 @echo "Generating MUI Fortran wrapper 3D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_f_wrapper_3d.cpp -o mui_f_wrapper_3d.o

mui_f_wrapper_general.mod: mui_f_wrapper_general.f90
	 @echo "Generating MUI Fortran wrapper general mod file..."
	 ${FC} ${FFLAGS} -c mui_f_wrapper_general.f90 -o mui_f_wrapper_general.mod

mui_f_wrapper_1d.mod: mui_f_wrapper_1d.f90
	 @echo "Generating MUI Fortran wrapper 1D mod file..."
	 ${FC} ${FFLAGS} -c mui_f_wrapper_1d.f90 -o mui_f_wrapper_1d.mod

mui_f_wrapper_2d.mod: mui_f_wrapper_2d.f90
	 @echo "Generating MUI Fortran wrapper 2D mod file..."
	 ${FC} ${FFLAGS} -c mui_f_wrapper_2d.f90 -o mui_f_wrapper_2d.mod

mui_f_wrapper_3d.mod: mui_f_wrapper_3d.f90
	 @echo "Generating MUI Fortran wrapper 3D mod file..."
	 ${FC} ${FFLAGS} -c mui_f_wrapper_3d.f90 -o mui_f_wrapper_3d.mod

test:
	 @echo "Compiling and linking MUI Fortran wrapper unit test code..."
	 ${FC} ${FFLAGS} mui_f_wrapper_1d.f90 mui_f_wrapper_2d.f90 mui_f_wrapper_3d.f90 unit_test.f90 -o unit_test_fortran_wrapper -L. -lstdc++ -lmpi_cxx -lmui_f_wrapper
	 ${FC} ${FFLAGS} mui_f_wrapper_1d.f90 mui_f_wrapper_2d.f90 mui_f_wrapper_3d.f90 unit_test_multi.f90 -o unit_test_multi_fortran_wrapper -L. -lstdc++ -lmpi_cxx -lmui_f_wrapper
	 @echo "Launching unit test code for direct interface creation..."
	 ${MPI} -np 1 ./unit_test_fortran_wrapper domain1 interface : -np 1 ./unit_test_fortran_wrapper domain2 interface
	 @echo "Launching unit test code for interface creation using helper fucntion..."
	 ${MPI} -np 1 ./unit_test_multi_fortran_wrapper domain1 interface 2 : -np 1 ./unit_test_multi_fortran_wrapper domain2 interface 2

clean:
	 rm -rf *.o *.mod *.so *.a unit_test_fortran_wrapper unit_test_multi_fortran_wrapper
