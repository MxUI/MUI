CC=mpicc
CCX=mpic++
MPI=mpirun

default: mui_3d.o
	ar rcs libwrapperC.a mui_3d.o

mui_3d.o: mui_3d.cpp unit_test.c
	@echo "Generating C-wrapper object file..."
	${CCX} -O3 -std=c++11 -c mui_3d.cpp -o mui_3d.o
	@echo "Compiling and linking unit test code..."
	${CC} -O3 -std=c11 unit_test.c mui_3d.o -o unit_test.x -lstdc++ -lmpi_cxx -lm
	@echo "Launching unit test code..."
	${MPI} -np 1 ./unit_test.x mpi://domain1/interface : -np 1 ./unit_test.x mpi://domain2/interface
	@echo "Usage: "
	@echo "Step 1) include mui_3d.h in C source"
	@echo "Step 2) Compile with a C compiler & link"

clean:
	rm -rf *.o *.a unit_test.x
