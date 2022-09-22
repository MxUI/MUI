# Define compilers and MPI run option
CC=mpicc
CCX=mpic++
MPI=mpirun
AR=ar

# C compilation flags
CFLAGS = -O3 -std=c11
# C++ compilation flags
CXXFLAGS = -O3 -std=c++11 -fpic

ifdef USE_RBF
    # Path to Eigen library
    EIGEN_PATH = /usr/include/eigen3
    # C compilation flags
    CFLAGS += -DUSE_RBF
    # C++ compilation flags
    CXXFLAGS += -DUSE_RBF -I${EIGEN_PATH}
endif

default: mui_c_wrapper_general.o mui_c_wrapper_1d.o mui_c_wrapper_2d.o mui_c_wrapper_3d.o
	 @echo "Compiling and linking MUI C wrapper..."
	 ${CCX} ${CFLAGS} -shared -o libmui_c_wrapper.so mui_c_wrapper_general.o mui_c_wrapper_1d.o mui_c_wrapper_2d.o mui_c_wrapper_3d.o
	 ${AR} rcs libmui_c_wrapper.a mui_c_wrapper_general.o mui_c_wrapper_1d.o mui_c_wrapper_2d.o mui_c_wrapper_3d.o

mui_c_wrapper_general.o: mui_c_wrapper_general.cpp
	 @echo "Generating C-wrapper general functions object file..."
	 ${CCX} ${CXXFLAGS} -c mui_c_wrapper_general.cpp -o mui_c_wrapper_general.o

mui_c_wrapper_1d.o: mui_c_wrapper_1d.cpp
	 @echo "Generating C-wrapper 1D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_c_wrapper_1d.cpp -o mui_c_wrapper_1d.o

mui_c_wrapper_2d.o: mui_c_wrapper_2d.cpp
	 @echo "Generating C-wrapper 2D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_c_wrapper_2d.cpp -o mui_c_wrapper_2d.o

mui_c_wrapper_3d.o: mui_c_wrapper_3d.cpp
	 @echo "Generating C-wrapper 3D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_c_wrapper_3d.cpp -o mui_c_wrapper_3d.o

test:
	 @echo "Compiling and linking MUI C wrapper unit test code..."
	 ${CC} ${CFLAGS} unit_test_single.c -o unit_test_c_wrapper_single -L${PWD} -lstdc++ -lmpi_cxx -lm -lmui_c_wrapper
	 ${CC} ${CFLAGS} unit_test_multi.c -o unit_test_c_wrapper_multi -L${PWD} -lstdc++ -lmpi_cxx -lm -lmui_c_wrapper
	 @echo "Launching C wrapper unit test code for single interface creation..."
	 ${MPI} -np 1 ./unit_test_c_wrapper_single domain1 interface : -np 1 ./unit_test_c_wrapper_single domain2 interface
	 @echo "Launching C wrapper unit test code for multi interface creation..."
	 ${MPI} -np 1 ./unit_test_c_wrapper_multi domain1 interface 2 : -np 1 ./unit_test_c_wrapper_multi domain2 interface 2

clean:
	 rm -rf *.o *.a *.so unit_test_c_wrapper_single unit_test_c_wrapper_multi
