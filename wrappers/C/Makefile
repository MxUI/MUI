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

default: mui_c_wrapper_1d.o
	 ${CCX} -shared -o lib_mui_c_wrapper.so mui_c_wrapper_1d.o
	 ${AR} rcs lib_mui_c_wrapper.a mui_c_wrapper_1d.o

mui_c_wrapper_1d.o: mui_c_wrapper_1d.cpp unit_test.c
	 @echo "Generating C-wrapper 1D object file..."
	 ${CCX} ${CXXFLAGS} -c mui_c_wrapper_1d.cpp -o mui_c_wrapper_1d.o
	 @echo "Compiling and linking unit test code..."
	 ${CC} ${CFLAGS} unit_test.c mui_c_wrapper_1d.o -o unit_test_c_wrapper -lstdc++ -lmpi_cxx -lm
	 @echo "Launching unit test code..."
	 ${MPI} -np 1 ./unit_test_c_wrapper mpi://domain1/interface : -np 1 ./unit_test_c_wrapper mpi://domain2/interface
	 @echo "Usage: "
	 @echo "Step 1) include mui_c_wapper_[1d/2d/3d].h in C source"
	 @echo "Step 2) Compile with a C compiler & link against library"

clean:
	 rm -rf *.o *.a *.so unit_test
