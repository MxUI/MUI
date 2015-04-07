default: mui_3df.cpp unit_test.f90
	@echo "Generating F-wrapper object file..."
	mpicxx -std=c++11 -c mui_3df.cpp -o mui_3df.o
	@echo "Compiling unit test code..."
	mpif90  -c unit_test.f90 -o unit_test.o
	mpic++ unit_test.o mui_3d.o -o unit_test.x
	@echo "Launching unit test code..."
	mpirun -np 1 ./unit_test.x mpi://domain1/interface : -np 1 ./unit_test.x mpi://domain2/interface
	@echo "Usage: "
	@echo "Step 1) Compile with a Fortran compiler but do not link"
	@echo "Step 2) Link with a C++ compiler the object code together with mui_3df.o"   
