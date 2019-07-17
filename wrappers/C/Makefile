default: mui_3d.o
	ar rcs libwrapperC.a mui_3d.o

mui_3d.o: mui_3d.cpp unit_test.c
	@echo "Generating C-wrapper object file..."
	mpic++ -std=c++11 -c mui_3d.cpp -o mui_3d.o
	@echo "Compiling and linking unit test code..."
	mpicc  unit_test.c mui_3d.o -o unit_test.x -lstdc++ -lmpi_cxx -lmpi -lm
	@echo "Launching unit test code..."
	mpirun -np 1 ./unit_test.x mpi://domain1/interface : -np 1 ./unit_test.x mpi://domain2/interface
	@echo "Usage: "
	@echo "Step 1) include mui_3d.h in C source"
	@echo "Step 2) Compile with a C compiler & link"
	
clean:
	rm *.o *.a unit_test.x