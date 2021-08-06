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
EIGEN_PATH = /home/stephen/eclipse-workspace/eigen-3.3.9
# C compilation flags
CFLAGS += -DUSE_RBF
# C++ compilation flags
CXXFLAGS += -DUSE_RBF -I${EIGEN_PATH}
endif

default: mui_c_wrapper_general.o mui_c_wrapper_1d.o mui_c_wrapper_2d.o mui_c_wrapper_3d.o
	 ${CCX} -shared -o libmui_c_wrapper.so mui_c_wrapper_1d.o mui_c_wrapper_2d.o mui_c_wrapper_3d.o
	 ${AR} rcs libmui_c_wrapper.a mui_c_wrapper_1d.o mui_c_wrapper_2d.o mui_c_wrapper_3d.o
	 @echo "Compiling and linking unit test code..."
	 ${CC} ${CFLAGS} unit_test.c -o unit_test_c_wrapper -L${PWD} -lstdc++ -lmpi_cxx -lm -lmui_c_wrapper
	 @echo "Launching unit test code..."
	 ${MPI} -np 1 ./unit_test_c_wrapper domain1 interface : -np 1 ./unit_test_c_wrapper domain2 interface
	 @echo "************************************************************"
	 @echo "*                    C wrapper usage                       *"
	 @echo "*                                                          *"
	 @echo "* Step 1) include mui_c_wapper_[1d/2d/3d].h in C source    *"
	 @echo "* Step 2) Compile with a C compiler & link against library *"
	 @echo "************************************************************"

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

clean:
	 rm -rf *.o *.a *.so unit_test_c_wrapper
