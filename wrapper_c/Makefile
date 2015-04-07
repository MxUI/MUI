default: mui_3d.cpp unit_test.c
	@echo "Generating C-wrapper object file..."
	mpic++ -std=c++11 -c mui_3d.cpp -o mui_3d.o
	@echo "Compiling unit test code..."
	mpicc  -c unit_test.c -o unit_test.o
	mpic++ unit_test.o mui_3d.o -o unit_test.x
	@echo "Launching unit test code..."
	mpirun -np 1 ./unit_test.x mpi://domain1/interface : -np 1 ./unit_test.x mpi://domain2/interface
	@echo "Usage: "
	@echo "Step 1) include mui_3d.h in C source"
	@echo "Step 2) Compile with a C compiler but do not link"
	@echo "Step 3) Link with a C++ compiler the object code together with mui_3d.o"   